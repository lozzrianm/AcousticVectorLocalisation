clc; clear all; close all;

% Single source localisation via ULA — experimental data (sequential recordings)
%
% Constructs an 8-element ULA from two sequential 4-mic recordings on a
% 6-channel data logger. The two recordings are windowed into snapshots
% independently, then the snapshot sets are pooled before forming the CSM.
% This avoids introducing an arbitrary inter-recording phase discontinuity
% into the coherence estimate.
%
% Off-diagonal CSM blocks (cross-terms between recordings) are zero by
% construction — the recordings are temporally independent so there is no
% coherent cross-correlation to estimate between them.
%
% Two independent calibration corrections are available, applied inside
% create_csm_pooled per frequency bin before CSM accumulation:
%   (1) GAIN CORRECTION — loaded from mic_gain_correction.mat
%       Linear gain multipliers per frequency bin, applied via interp1.
%   (2) PHASE CORRECTION — defined in the calibration section below
%       Phase offsets (rad) measured at ref_freq, converted to time delays
%       and applied as exp(-1i * 2*pi*f * tau) at each bin.
%
% In addition, an optional fixed-gain MEMS calibration block runs
% before CSM construction, applying linear per-channel multipliers
% derived from pairwise pistonphone measurements. This is independent
% of the frequency-table gain correction above.
%
% Written by L Marshall 30/03/2026
% Frequency-locked bin selection, grid override, MEMS gain block
% by L Marshall 29/04/2026
%
% Changes (29/04/2026):
%   - FFT length and bin selection now centred on the source frequency
%     rather than driven by the spatial-resolution lower-bound formula.
%     Analysis bandwidth is a narrow window around f_test where the
%     signal energy actually lives — appropriate for performance
%     characterisation of a known monochromatic test source. The previous
%     broadband scheme integrated many noise-only bins, biasing results.
%   - grid_resolution_fixed parameter added — set to a positive value to
%     override the frequency-dependent λ/N grid with a constant step,
%     useful for isolating grid quantisation effects from genuine
%     frequency-dependent performance variation.
%   - MEMS per-channel amplitude calibration block added in the signal
%     conditioning section, parallel to the AVA experimental script.


%% DEFINE INPUT VARIABLES %%

% Array 1 input file parameters
wav_filename_a1 = '20260330-015356(UTC)-Trial12.MEMS.ULA.4MIC.0.4m.5.143cmdelta.pos1fixed.1khz-0166476623.wav';
csv_filename_a1 = 'MPA416_MEMS_Array1_Signal_fixed.csv';
convert_wav_a1 = true;

% Array 2 input file parameters
wav_filename_a2 = '20260330-014125(UTC)-Trial11.MEMS.ULA.4MIC.0.4m.5.143cmdelta.pos2fixed.1khz-0165727396.wav';
csv_filename_a2 = 'MPA416_MEMS_Array2_Signal_fixed.csv';
convert_wav_a2 = true;

% Physical environment parameters
c_0 = 340; %speed of sound (m/s)

% Expected source position and frequency — for reference and error metrics
source_x = -0.4;
source_y = 0;
source_positions = [source_x, source_y, 0];
source_frequencies = [1000]; %source frequency (Hz)

% Array geometry
N_mics_per_recording = 4; %microphones per recording
N_mics = 8; %total microphones in combined ULA
x_a = 0; %array x coordinate (m)
z_a = 0; %array z coordinate (m)

d_y = 0.36 / 7; %inter-element spacing (m) — 0.051429 m
y_centre_a1 = 0.02 - 1.5 * d_y; %sub-array 1 centre (m) — -0.057143 m
y_centre_a2 = -0.34 + 1.5 * d_y; %sub-array 2 centre (m) — -0.262857 m
y_a = (y_centre_a1 + y_centre_a2) / 2; %combined geometric centre (m) — -0.160000 m

% Processing parameters
overlap = 0.5; %window overlap fraction
loading = 1e-4; %regularisation parameter
analysis_halfwidth_hz = 10; %half-width of analysis band around source frequency (Hz)
target_binwidth_hz = 1; %target FFT bin spacing within analysis band (Hz)

% Grid parameters — frequency-dependent by default, with a fixed override
% available to isolate grid quantisation effects
grid_pts_per_lambda = 20; %grid points per wavelength
grid_resolution_fixed = []; %override grid step (m); [] = use λ/grid_pts_per_lambda
margin_fixed = 0.5; %fixed physical search margin (m)


%% CALIBRATION CONFIGURATION %%

% Toggle flags
apply_gain_correction = false;
apply_phase_correction = false;

% Gain correction mat file
gain_corr_mat_file = 'mems_mic_gain_correction.mat';

% Phase correction — offsets in radians at ref_freq
% Same physical microphones used for both recordings
% Sign convention: positive = channel leads ch 1, correction lags
% Time-delay model: phi_at_f = phi_at_ref * (f / ref_freq), valid when the
% physical source of the offset is frequency-flat (cable/ADC group delay)
ref_freq = 1000;

phase_offset_rad = [
    +0.000000; %ch 1
    +0.216973; %ch 2
    +0.039670; %ch 3
    +0.196226; %ch 4
];

% MEMS per-channel amplitude calibration — applied directly to the time
% series before snapshot construction. Multipliers derived from pairwise
% pistonphone Vpp measurements at 1 kHz:
%   M2 vs M1: 312/324 = 0.9630  -> correction 1.0385
%   M3 vs M1: 308/324 = 0.9506  -> correction 1.0519
%   M4 vs M1: 279/306 = 0.9118  -> correction 1.0968
% Same physical microphones are used for both recordings, so the same
% multipliers apply to both Array 1 and Array 2 time series.
apply_mems_gain_cal = false;
mems_gain_correction = [1.0000; 1.0385; 1.0519; 1.0968];


%% LOAD GAIN CORRECTION TABLE %%

if apply_gain_correction
    if exist(gain_corr_mat_file, 'file')
        load(gain_corr_mat_file, 'gain_corr_table', 'corr_freq_vec', 'phase_offset_rad');
        fprintf('Loaded gain correction table: [%d freqs x %d channels]\n', ...
            size(gain_corr_table, 1), size(gain_corr_table, 2));
        fprintf('  Frequency range: %.0f - %.0f Hz\n', ...
            min(corr_freq_vec), max(corr_freq_vec));
        fprintf('Loaded phase_offset_rad from mat file: ');
        fprintf('Ch%d = %+.4f  ', [(1:length(phase_offset_rad)); phase_offset_rad(:)']);
        fprintf('\n');
    else
        warning('mic_gain_correction.mat not found. Gain correction disabled.');
        apply_gain_correction = false;
        gain_corr_table = [];
        corr_freq_vec = [];
    end
else
    gain_corr_table = [];
    corr_freq_vec = [];
end

fprintf('\n<strong>Calibration Configuration</strong>\n');
fprintf('Gain correction (table): %s\n', string(apply_gain_correction));
fprintf('Phase correction:        %s\n', string(apply_phase_correction));
fprintf('MEMS fixed gain cal:     %s\n', string(apply_mems_gain_cal));
if apply_phase_correction
    fprintf('Phase reference frequency: %d Hz\n', ref_freq);
    fprintf('Phase offsets (rad): ');
    fprintf('Ch%d = %+.4f  ', [(1:length(phase_offset_rad)); phase_offset_rad(:)']);
    fprintf('\n');
    tau_vec = phase_offset_rad / (2 * pi * ref_freq);
    fprintf('Implied time delays (us): ');
    fprintf('Ch%d = %+.2f  ', [(1:length(tau_vec)); tau_vec(:)' * 1e6]);
    fprintf('\n');
end
fprintf('\n');


%% CONVERT WAV FILES TO CSV %%

% Array 1
if convert_wav_a1
    fprintf('<strong>Converting Array 1 WAV to CSV</strong>\n');
    if ~exist(wav_filename_a1, 'file')
        error('Array 1 WAV file not found: %s', wav_filename_a1);
    end
    [audio_a1, fs_a1] = audioread(wav_filename_a1);
    fprintf('  Array 1: %d samples, %d channels, %d Hz\n', ...
        size(audio_a1, 1), size(audio_a1, 2), fs_a1);
    time_a1 = (0:size(audio_a1, 1) - 1)' / fs_a1;
    writematrix([time_a1, audio_a1], csv_filename_a1);
    fprintf('  Saved: %s\n\n', csv_filename_a1);
end

% Array 2
if convert_wav_a2
    fprintf('<strong>Converting Array 2 WAV to CSV</strong>\n');
    if ~exist(wav_filename_a2, 'file')
        error('Array 2 WAV file not found: %s', wav_filename_a2);
    end
    [audio_a2, fs_a2] = audioread(wav_filename_a2);
    fprintf('  Array 2: %d samples, %d channels, %d Hz\n', ...
        size(audio_a2, 1), size(audio_a2, 2), fs_a2);
    time_a2 = (0:size(audio_a2, 1) - 1)' / fs_a2;
    writematrix([time_a2, audio_a2], csv_filename_a2);
    fprintf('  Saved: %s\n\n', csv_filename_a2);
end


%% LOAD SIGNAL DATA FROM CSV %%

fprintf('<strong>Loading signal data</strong>\n');

% Array 1
if ~exist(csv_filename_a1, 'file')
    error('Array 1 CSV not found: %s\nRun WAV conversion first.', csv_filename_a1);
end
data_a1 = readmatrix(csv_filename_a1);
time_a1 = data_a1(:, 1);
tx_a1 = data_a1(:, 2:end).'; %rows = sensors, columns = time samples

if size(tx_a1, 1) ~= N_mics_per_recording
    warning('Array 1: expected %d microphones but found %d in CSV.', ...
        N_mics_per_recording, size(tx_a1, 1));
end
fprintf('  Array 1: %d samples from %d microphones (%.3f s)\n', ...
    size(tx_a1, 2), size(tx_a1, 1), time_a1(end));

% Array 2
if ~exist(csv_filename_a2, 'file')
    error('Array 2 CSV not found: %s\nRun WAV conversion first.', csv_filename_a2);
end
data_a2 = readmatrix(csv_filename_a2);
time_a2 = data_a2(:, 1);
tx_a2 = data_a2(:, 2:end).';

if size(tx_a2, 1) ~= N_mics_per_recording
    warning('Array 2: expected %d microphones but found %d in CSV.', ...
        N_mics_per_recording, size(tx_a2, 1));
end
fprintf('  Array 2: %d samples from %d microphones (%.3f s)\n\n', ...
    size(tx_a2, 2), size(tx_a2, 1), time_a2(end));


%% SIGNAL CONDITIONING %%

tx_a1 = tx_a1 - mean(tx_a1, 2); %remove DC offset
tx_a2 = tx_a2 - mean(tx_a2, 2);

% TEMPORARY DIAGNOSTIC — ULA channel reorder for end-to-end flip + middle swap
% Diagnostic FFT@1kHz showed CSV channel 4 sits at the most positive y position,
% then ch2, ch3, ch1 in descending y. Apply the same mapping to both recordings
% (same physical wiring used for both).
% So in slot order [row1, row2, row3, row4] (descending y) we want CSV
% columns [ch4, ch2, ch3, ch1].

% UNCOMMENT THESE IF USING THE MEMS TRIALS FROM MARCH 

tx_a1 = tx_a1([4 2 3 1], :);
tx_a2 = tx_a2([4 2 3 1], :);

% MEMS PER-CHANNEL AMPLITUDE CALIBRATION %
% Applied as a fixed linear gain on the time series before snapshot
% construction. Distinct from the frequency-dependent gain table loaded
% from mems_mic_gain_correction.mat — the two paths are independently
% toggled and operate at different stages of the pipeline.
if apply_mems_gain_cal
    tx_a1 = tx_a1 .* mems_gain_correction;
    tx_a2 = tx_a2 .* mems_gain_correction;
    fprintf('  Applied MEMS gain calibration (M2: +%.3f dB, M3: +%.3f dB, M4: +%.3f dB)\n\n', ...
        20*log10(mems_gain_correction(2)), 20*log10(mems_gain_correction(3)), ...
        20*log10(mems_gain_correction(4)));
end

% Channel order preserved as recorded — mic 1 is the most positive y element
% within each sub-array, consistent with the physical labelling convention


%% ARRAY GEOMETRY SETUP %%

fprintf('<strong>Setting up combined ULA geometry</strong>\n');

% Descending y order — row 1 of tx_a1/tx_a2 is the most positive y element
n_vec = linspace((N_mics_per_recording-1)/2, -(N_mics_per_recording-1)/2, N_mics_per_recording);
y_mics_a1 = y_centre_a1 + d_y * n_vec;
y_mics_a2 = y_centre_a2 + d_y * n_vec;

% Combined mic positions: array 1 (upper y, rows 1-4) then array 2 (lower y, rows 5-8)
% Row order matches the signal row order after concatenation and CSM block assignment
y_mics_combined = [y_mics_a1, y_mics_a2];
mic_positions = [x_a * ones(1, N_mics); y_mics_combined; z_a * ones(1, N_mics)];

fprintf('  Array 1 mic positions (y): ');
fprintf('%.4f  ', y_mics_a1);
fprintf('m\n');
fprintf('  Array 2 mic positions (y): ');
fprintf('%.4f  ', y_mics_a2);
fprintf('m\n');
fprintf('  Combined ULA (%d mics, y): ', N_mics);
fprintf('%.4f  ', y_mics_combined);
fprintf('m\n');
fprintf('  Total aperture: %.4f m\n', max(y_mics_combined) - min(y_mics_combined));
fprintf('  Geometric centre: y = %.4f m\n\n', mean(y_mics_combined));


%% FREQUENCY DOMAIN SETUP %%

% Both recordings assumed at the same sample rate — use Array 1 time vector
d_t = time_a1(2) - time_a1(1);
F_s = 1 / d_t;

fprintf('<strong>Frequency domain setup</strong>\n');
fprintf('  Sample rate: %.0f Hz | Time step: %.6f s\n', F_s, d_t);


%% CHOOSE FFT LENGTH AND ANALYSIS BINS %%
% Pick size_fft so that the source frequency lands on a bin centre at
% the target spacing — keeps signal energy in a small number of adjacent
% bins rather than smearing it across the spectrum, and avoids integrating
% noise-only bins into the broadband MVDR sum.

f_test = source_frequencies(1);
f_band_lo = f_test - analysis_halfwidth_hz;
f_band_hi = f_test + analysis_halfwidth_hz;

cycles_per_window = round(f_test / target_binwidth_hz);
size_fft = round(cycles_per_window * F_s / f_test);
actual_binwidth = F_s / size_fft;

% Cap at recording length to maintain reasonable snapshot count
N_a1 = size(tx_a1, 2);
N_a2 = size(tx_a2, 2);
N_min = min(N_a1, N_a2);
if size_fft > N_min / 4
    size_fft = 2^floor(log2(N_min / 4));
    actual_binwidth = F_s / size_fft;
    fprintf('  Capped size_fft at %d (recording too short for target binwidth)\n', size_fft);
end

fft_vec = F_s * (0:(size_fft - 1)) / size_fft;

bin_mask = (fft_vec >= f_band_lo) & (fft_vec <= f_band_hi);
bin_index = find(bin_mask);
bin_freqs = fft_vec(bin_index);
num_bins = length(bin_index);

if num_bins == 0
    [~, nearest] = min(abs(fft_vec(1:floor(size_fft/2)) - f_test));
    bin_index = nearest;
    bin_freqs = fft_vec(nearest);
    num_bins = 1;
    fprintf('  Analysis band narrower than bin width — using single nearest bin\n');
end

fprintf('  size_fft = %d | bin width = %.3f Hz | %d bins (%.1f-%.1f Hz)\n\n', ...
    size_fft, actual_binwidth, num_bins, min(bin_freqs), max(bin_freqs));


%% PLOT FREQUENCY SPECTRA %%

figure('Name', 'Frequency Spectra - Recording 1 (Array 1 position)');
f_plot = F_s * (0:floor(N_a1/2)) / N_a1;
hold on;
for sensor = 1:N_mics_per_recording
    Y = fft(tx_a1(sensor, :));
    Y = Y(1:floor(N_a1/2) + 1);
    Y_mag = abs(Y) / N_a1;
    Y_mag(2:end-1) = 2 * Y_mag(2:end-1);
    plot(f_plot, 20 * log10(Y_mag + eps), 'DisplayName', sprintf('Sensor %d', sensor));
end
hold off;
xlim([0 max(2 * f_test, 2000)]);
xlabel('Frequency (Hz)', 'FontName', 'Times New Roman', 'FontSize', 14);
ylabel('Amplitude (dB)', 'FontName', 'Times New Roman', 'FontSize', 14);
legend('show', 'Location', 'best');
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);

figure('Name', 'Frequency Spectra - Recording 2 (Array 2 position)');
f_plot = F_s * (0:floor(N_a2/2)) / N_a2;
hold on;
for sensor = 1:N_mics_per_recording
    Y = fft(tx_a2(sensor, :));
    Y = Y(1:floor(N_a2/2) + 1);
    Y_mag = abs(Y) / N_a2;
    Y_mag(2:end-1) = 2 * Y_mag(2:end-1);
    plot(f_plot, 20 * log10(Y_mag + eps), 'DisplayName', sprintf('Sensor %d', sensor));
end
hold off;
xlim([0 max(2 * f_test, 2000)]);
xlabel('Frequency (Hz)', 'FontName', 'Times New Roman', 'FontSize', 14);
ylabel('Amplitude (dB)', 'FontName', 'Times New Roman', 'FontSize', 14);
legend('show', 'Location', 'best');
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);


%% CREATE SNAPSHOTS AND POOLED CROSS-SPECTRAL MATRIX %%

% Snapshots created independently from each recording to avoid inter-recording
% phase discontinuity. The two snapshot sets are pooled — equivalent to
% averaging the per-recording CSMs weighted by snapshot count.
%
% The resulting N_mics x N_mics CSM has block structure:
%   rows/cols 1:N_per     = Array 1 position mics (upper y, most positive y first)
%   rows/cols N_per+1:end = Array 2 position mics (lower y, most positive y first)
%
% Off-diagonal blocks are zero — recordings are temporally independent.

fprintf('<strong>Creating snapshots and pooled cross-spectral matrix</strong>\n');

window = hanning(size_fft)';

snapshots_a1 = make_snapshots(tx_a1, size_fft, overlap, window);
snapshots_a2 = make_snapshots(tx_a2, size_fft, overlap, window);

fprintf('  Recording 1: %d snapshots\n', size(snapshots_a1, 3));
fprintf('  Recording 2: %d snapshots\n', size(snapshots_a2, 3));

r = create_csm_pooled(snapshots_a1, snapshots_a2, bin_index, N_mics, N_mics_per_recording, ...
    num_bins, bin_freqs, apply_gain_correction, gain_corr_table, corr_freq_vec, ...
    apply_phase_correction, phase_offset_rad, ref_freq);

fprintf('  Pooled CSM: %dx%dx%d\n\n', size(r));

%% (after r is built and mic_positions is set up — i.e. after STEP 8 / before STEP 9) %%

% SUB-ARRAY SELECTION FOR 4-MIC ULA BASELINE COMPARISON %
% Selects 4 of the 8 mics for a reduced-aperture ULA test. Pattern A is
% an alternating uniform ULA at 0.103 m spacing, full-bay aperture.
% Set use_subarray = false to keep the full 8-mic array.
use_subarray = true;
subarray_pattern = 'A';  % 'A' = uniform 4-mic alternating; 'C' = sparse non-uniform full-aperture

if use_subarray
    y_a = mean(mic_positions(2,:));
    fprintf('  Updated array centre y: %.4f m\n', y_a);
    switch subarray_pattern
        case 'A'  % rows 1, 3, 5, 7 — uniform 0.103 m spacing across 0.309 m aperture
            sel = [1 3 5 7];
            label = 'Pattern A (uniform 4-mic alternating, 0.309 m aperture)';
        case 'C'  % rows 1, 3, 6, 8 — non-uniform sparse, 0.360 m aperture
            sel = [1 3 6 8];
            label = 'Pattern C (sparse non-uniform 4-mic, 0.360 m aperture)';
        otherwise
            error('Unknown subarray_pattern: %s', subarray_pattern);
    end

    fprintf('\n<strong>Sub-array selection: %s</strong>\n', label);
    fprintf('  Selected rows: ');
    fprintf('%d ', sel);
    fprintf('\n  Selected y positions: ');
    fprintf('%+.4f ', mic_positions(2, sel));
    fprintf('m\n');

    r = r(sel, sel, :);
    mic_positions = mic_positions(:, sel);
    N_mics = length(sel);
    fprintf('  Reduced CSM: %dx%dx%d\n\n', size(r));
end
%% DEFINE GRID SEARCH AREA %%

fprintf('<strong>Defining search area</strong>\n');

system_centre_y = -0.16; %geometric midpoint of full aperture (m)
system_centre_y = y_a; %geometric midpoint of (possibly reduced) aperture (m)
[x_scan, y_scan, X_grid, Y_grid, candidate_points, grid_res] = ...
    build_search_grid(f_test, c_0, 0.4, ...
    system_centre_y, grid_pts_per_lambda, margin_fixed, grid_resolution_fixed);

fprintf('  Search grid: %d x %d points, resolution: %.4f m\n\n', ...
    length(x_scan), length(y_scan), grid_res);


%% PERFORM MVDR BEAMFORMING %%

fprintf('<strong>Running MVDR Beamforming</strong>\n');
fprintf('  Regularisation parameter: %.2e\n\n', loading);

response_db = mvdr_beamforming(r, mic_positions, candidate_points, bin_freqs, c_0, loading);


%% 1D SCAN %%

fprintf('1D Scan (x = %.3f m):\n', source_x);
plot_1dscan_mvdr(r, mic_positions, source_x, source_y, bin_freqs, c_0, loading);


%% 2D SCAN %%

fprintf('2D Scan:\n');
[est_x, est_y] = plot_2dscan_mvdr(r, mic_positions, candidate_points, bin_freqs, c_0, ...
    y_scan, x_scan, source_x, source_y, X_grid, Y_grid, loading, grid_res);


%% PERFORMANCE METRICS %%

fprintf('\n<strong>Performance Metrics:</strong>\n');

est_pos = [est_x, est_y];
true_pos = [source_x, source_y];
errors = true_pos - est_pos;

radial_error = norm(errors);
MSE = sum(errors.^2);
MSE_percent = (MSE / (grid_res^2)) * 100;

true_angle = atan2(source_y - y_a, source_x - x_a);
est_angle = atan2(est_y - y_a, est_x - x_a);
angular_error_deg = rad2deg(abs(true_angle - est_angle));
if angular_error_deg > 180
    angular_error_deg = 360 - angular_error_deg;
end

fprintf('  Estimated position: (%.4f, %.4f) m\n', est_x, est_y);
fprintf('  Position error: dx = %.4f m, dy = %.4f m\n', errors(1), errors(2));
fprintf('  Radial error: %.4f m\n', radial_error);
fprintf('  Mean squared error: %.4f m^2\n', MSE);
fprintf('  MSE (normalised): %.3f%% of grid resolution\n', MSE_percent);
fprintf('  Angular error: %.4f deg\n', angular_error_deg);


%% BEAM PATTERN ANALYSIS %%

fprintf('\n<strong>GENERATING BEAM PATTERNS</strong>\n');

array_centre_2d = [x_a; y_a];
distances = sqrt(sum((source_positions(:,1:2) - array_centre_2d').^2, 2));
radius = median(distances);
fprintf('Using radius = %.3f m (median source distance)\n', radius);

[~, ~, ~] = beam_pattern_ula(...
    r, mic_positions, array_centre_2d, source_positions(:,1:2), ...
    bin_freqs, c_0, loading, radius, ...
    sprintf('%d-Mic ULA Beam Pattern (MVDR)', N_mics));


%% FUNCTION DEFINITIONS %%


% FUNCTION: MAKE_SNAPSHOTS
% Make time-domain snapshots of input signal
function snapshots = make_snapshots(tx, size_fft, overlap, window)

    step = round(overlap * size_fft);
    num_snap = floor((size(tx, 2) - size_fft) / step) + 1;

    start_idx = 1 + (0:(num_snap - 1)) * step;
    idx_matrix = start_idx + (0:size_fft - 1)';

    snapshots = tx(:, idx_matrix);
    snapshots = reshape(snapshots, size(tx, 1), size_fft, num_snap);
    snapshots = snapshots .* reshape(window, [1, size_fft, 1]);
end


% FUNCTION: CREATE_CSM_POOLED
% Pooled cross-spectral matrix from two independent recordings
%
% Snapshots from each recording are transformed to the frequency domain
% independently. Gain and phase corrections are applied per recording before
% the outer product is formed. The two partial CSMs are then averaged,
% weighted by their respective snapshot counts.
%
% The resulting N_mics x N_mics CSM has block structure:
%   rows/cols 1:N_per     = Array 1 position mics (upper y, most positive y first)
%   rows/cols N_per+1:end = Array 2 position mics (lower y, most positive y first)
%
% Off-diagonal blocks are zero — recordings are temporally independent.
%
% Gain correction: gain_corr_table interpolated at each bin frequency and
% applied column-wise to frequency-domain snapshots.
%
% Phase correction: phase_offset_rad (measured at ref_freq) converted to
% time delay tau = phi / (2*pi*ref_freq), applied as exp(-1i*2*pi*f*tau).
function r = create_csm_pooled(snapshots_a1, snapshots_a2, bin_index, N_mics, N_per, ...
    num_bins, bin_freqs, apply_gain, gain_corr_table, corr_freq_vec, ...
    apply_phase, phase_offset_rad, ref_freq)

    r = zeros(N_mics, N_mics, num_bins);
    num_snap1 = size(snapshots_a1, 3);
    num_snap2 = size(snapshots_a2, 3);

    fx_a1 = fft(snapshots_a1, [], 2);
    fx_a2 = fft(snapshots_a2, [], 2);
    fx_a1 = fx_a1(:, bin_index, :);
    fx_a2 = fx_a2(:, bin_index, :);

    for jf = 1:num_bins
        freq = bin_freqs(jf);

        fb1 = squeeze(fx_a1(:, jf, :)); %(N_per x num_snap)
        fb2 = squeeze(fx_a2(:, jf, :));

        % Gain correction — applied independently to each recording
        if apply_gain && ~isempty(gain_corr_table)
            gain_vec = interp1(corr_freq_vec, gain_corr_table, freq, 'linear', 'extrap');
            fb1 = fb1 .* gain_vec(:);
            fb2 = fb2 .* gain_vec(:);
        end

        % Phase correction — applied independently to each recording
        if apply_phase && ~isempty(phase_offset_rad)
            tau_vec = phase_offset_rad(:) / (2 * pi * ref_freq);
            phase_corr_vec = exp(-1i * 2 * pi * freq * tau_vec);
            fb1 = fb1 .* phase_corr_vec;
            fb2 = fb2 .* phase_corr_vec;
        end

        R1 = (fb1 * fb1') / num_snap1;
        R2 = (fb2 * fb2') / num_snap2;

        % Block-diagonal pooled CSM
        r_bin = zeros(N_mics, N_mics);
        r_bin(1:N_per, 1:N_per) = R1;
        r_bin(N_per+1:N_mics, N_per+1:N_mics) = R2;

        r(:,:,jf) = r_bin;
    end
end


% FUNCTION: MVDR_BEAMFORMING
% MVDR beamforming with adaptive CBF-based diagonal loading
function response_db = mvdr_beamforming(r, mic_positions, candidate_points, bin_freqs, c_0, loading)

    num_points = size(candidate_points, 1);
    num_bins = size(r, 3);
    mvdr_responses = zeros(num_points, num_bins);

    for jf = 1:num_bins
        k = 2 * pi * bin_freqs(jf) / c_0;
        rf = squeeze(r(:, :, jf));

        [ur, er] = eig(rf);
        erv = real(diag(er));
        erv = max(erv, 1e-12);

        for n = 1:num_points
            l_pt = candidate_points(n, :).';
            r_lm = sqrt(sum((mic_positions - l_pt).^2, 1));
            v1 = exp(-1i * k * r_lm).';

            cbf_out = real(v1' * rf * v1);
            lambda = loading * cbf_out;
            if lambda == 0
                lambda = loading;
            end

            rxv = (ur * diag(1 ./ (erv + lambda)) * ur') * v1;
            denominator = v1' * rxv;

            if abs(denominator) < 1e-12
                mvdr_responses(n, jf) = 0;
            else
                vmvdr = rxv / denominator;
                mvdr_responses(n, jf) = abs(vmvdr' * rf * vmvdr);
            end
        end
    end

    responses_sum = sum(mvdr_responses, 2);
    response_db = 10 * log10(responses_sum + eps);
    response_db = response_db - max(response_db);
end


% FUNCTION: PLOT_2DSCAN_MVDR
% 2D MVDR beamformer scan — publication style
% Returns estimated source position from peak response
function [est_x, est_y] = plot_2dscan_mvdr(r, mic_positions, candidate_points, ...
    bin_freqs, c_0, y_scan, x_scan, source_x, source_y, X_grid, Y_grid, loading, grid_res)

    response_db = mvdr_beamforming(r, mic_positions, candidate_points, bin_freqs, c_0, loading);
    grid_response = reshape(response_db, length(y_scan), length(x_scan));

    figure('Color', 'w', 'Position', [100 100 520 440]);
    imagesc(x_scan, y_scan, grid_response);
    axis xy; axis tight;
    set(gca, 'PositionConstraint', 'innerposition');
    colormap(gca, 'jet');

    cb = colorbar;
    cb.FontName = 'Times New Roman';
    cb.FontSize = 11;
    cb.Label.String = 'Normalised Response (dB)';
    cb.Label.FontName = 'Times New Roman';
    cb.Label.FontSize = 12;
    cb.Label.Interpreter = 'latex';

    xlabel('$x$ (m)', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('$y$ (m)', 'Interpreter', 'latex', 'FontSize', 14);
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 12, ...
        'LineWidth', 1.0, 'Box', 'on', 'Layer', 'top');
    grid off; hold on;

    [est_x, est_y, ~] = refine_peak_2d(grid_response, x_scan, y_scan);

    fprintf('  Max response at (%.3f, %.3f) m\n', est_x, est_y);

    % True source — black cross
    plot(source_x, source_y, 'kx', ...
        'MarkerSize', 7, 'LineWidth', 1.5, 'DisplayName', 'True source');

    % Estimated source — black hollow circle
    plot(est_x, est_y, 'ko', ...
        'MarkerSize', 10, 'LineWidth', 1.5, 'MarkerFaceColor', 'none', ...
        'DisplayName', 'Estimated source');

    % Array aperture — vertical line with endpoint ticks
    y_aperture_top = max(mic_positions(2,:));
    y_aperture_bot = min(mic_positions(2,:));
    x_array = mean(mic_positions(1,:));
    tick_half = 0.05; %half-width of endpoint ticks (m)
    col_array = [0.7 0.7 0.7]; %medium grey

    plot([x_array, x_array], [y_aperture_bot, y_aperture_top], '-', ...
        'LineWidth', 1.5, 'Color', col_array, 'DisplayName', 'Array');
    plot([x_array - tick_half, x_array + tick_half], [y_aperture_top, y_aperture_top], '-', ...
        'LineWidth', 1.5, 'Color', col_array, 'HandleVisibility', 'off');
    plot([x_array - tick_half, x_array + tick_half], [y_aperture_bot, y_aperture_bot], '-', ...
        'LineWidth', 1.5, 'Color', col_array, 'HandleVisibility', 'off');

    legend('Location', 'northeast', 'FontName', 'Times New Roman', ...
        'FontSize', 10, 'Color', 'w', 'TextColor', 'k', 'Box', 'on');
end


% FUNCTION: PLOT_1DSCAN_MVDR
% 1D MVDR scan along y-axis at fixed x equal to source x coordinate
function plot_1dscan_mvdr(r, mic_positions, source_x, source_y, bin_freqs, c_0, loading)

    y_scan_line = linspace(source_y - 0.5, source_y + 0.5, 200);
    candidate_points = [source_x * ones(length(y_scan_line), 1), ...
                        y_scan_line(:), zeros(length(y_scan_line), 1)];

    response_line = mvdr_beamforming(r, mic_positions, candidate_points, bin_freqs, c_0, loading);

    [max_response, max_idx] = max(response_line);
    y_est = y_scan_line(max_idx);

    fprintf('  Max response: %.2f dB at y = %.3f m\n', max_response, y_est);

    figure('Name', '1D MVDR Scan');
    plot(y_scan_line, response_line, 'LineWidth', 2);
    xlabel('$y$ position (m)', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('Beamformer response (dB)', 'FontName', 'Times New Roman', 'FontSize', 14);
    grid on; hold on;
    xline(source_y, 'r--', 'LineWidth', 1.5, 'DisplayName', 'True source');
    legend('Beamformer response', 'True source', 'FontName', 'Times New Roman');
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
end


% FUNCTION: BEAM_PATTERN_ULA
% Near-field beam pattern for combined ULA MVDR
% Evaluates MVDR response at each azimuth at fixed range radius
function [theta_deg, beam_pattern, fig_handle] = beam_pattern_ula(...
    r, mic_positions, array_centre, source_positions, bin_freqs, c_0, ...
    loading, radius, plot_title)

    fprintf('\n<strong>Computing Near-Field Beam Pattern (ULA MVDR):</strong>\n');
    fprintf('  Range: %.3f m\n', radius);
    fprintf('  Array centre: (%.3f, %.3f) m\n', array_centre);

    num_angles = 360;
    theta_deg = linspace(0, 360, num_angles + 1);
    theta_deg = theta_deg(1:end-1); %remove duplicate at 360 deg
    theta_rad = deg2rad(theta_deg);
    num_bins = size(r, 3);
    mvdr_responses = zeros(num_angles, num_bins);

    for jf = 1:num_bins
        k = 2 * pi * bin_freqs(jf) / c_0;
        rf = squeeze(r(:, :, jf));

        [ur, er] = eig(rf);
        erv = real(diag(er));
        erv = max(erv, 1e-12);

        for angle_idx = 1:num_angles
            % 0 deg is along +Y axis (north), increases clockwise
            x_pos = array_centre(1) + radius * sin(theta_rad(angle_idx));
            y_pos = array_centre(2) + radius * cos(theta_rad(angle_idx));
            source_pos = [x_pos; y_pos; 0];

            r_lm = sqrt(sum((mic_positions - source_pos).^2, 1));
            v1 = exp(-1i * k * r_lm).';

            cbf_out = real(v1' * rf * v1);
            lambda = loading * cbf_out;
            if lambda == 0
                lambda = loading;
            end

            rxv = (ur * diag(1 ./ (erv + lambda)) * ur') * v1;
            denominator = v1' * rxv;

            if abs(denominator) > 1e-12
                vmvdr = rxv / denominator;
                mvdr_responses(angle_idx, jf) = abs(vmvdr' * rf * vmvdr);
            else
                mvdr_responses(angle_idx, jf) = 0;
            end
        end

        if mod(jf, 10) == 0 || jf == num_bins
            fprintf('  Processed %d/%d frequency bins\r', jf, num_bins);
        end
    end
    fprintf('\n');

    beam_pattern_linear = sum(mvdr_responses, 2);
    beam_pattern = 10 * log10(beam_pattern_linear + eps);
    beam_pattern = beam_pattern - max(beam_pattern);

    fig_handle = figure('Name', plot_title, 'Position', [100, 100, 700, 700]);
    polarplot(theta_rad, beam_pattern, 'b-', 'LineWidth', 2);
    ax = gca;
    ax.ThetaDir = 'clockwise';
    ax.ThetaZeroLocation = 'top'; %0 deg is North (+Y direction)
    rlim([min(beam_pattern), 0]);

    rticks_vals = linspace(min(beam_pattern), 0, 5);
    rticks(rticks_vals);
    rticklabels(arrayfun(@(x) sprintf('%.1f dB', x), rticks_vals, 'UniformOutput', false));

    if ~isempty(source_positions)
        hold on;
        for src = 1:size(source_positions, 1)
            source_vec = source_positions(src,:)' - array_centre;
            source_angle_rad = atan2(source_vec(1), source_vec(2));
            r_lim = rlim;
            polarplot([source_angle_rad, source_angle_rad], r_lim, '--r', 'LineWidth', 2, ...
                'DisplayName', sprintf('Source %d', src));
        end
        legend('Location', 'northoutside', 'Orientation', 'horizontal');
        hold off;
    end

    title(sprintf('%s (r = %.3f m)', plot_title, radius), ...
        'FontName', 'Times New Roman', 'FontSize', 16, 'FontWeight', 'bold');
    ax.FontName = 'Times New Roman';
    ax.FontSize = 12;

    fprintf('\n<strong>Beam Pattern Statistics:</strong>\n');
    [~, peak_idx] = max(beam_pattern);
    fprintf('  Peak response: 0.00 dB at %.1f degrees\n', theta_deg(peak_idx));

    threshold_3db = max(beam_pattern) - 3;
    above = beam_pattern(:) >= threshold_3db;
    d_above = diff([0; above; 0]);
    starts = find(d_above == 1);
    ends_arr = find(d_above == -1) - 1;

    main_lobe_region = find(starts <= peak_idx & ends_arr >= peak_idx, 1);
    if ~isempty(main_lobe_region)
        bw = theta_deg(ends_arr(main_lobe_region)) - theta_deg(starts(main_lobe_region));
        if bw < 0; bw = bw + 360; end
        fprintf('  -3 dB beamwidth: %.1f degrees (%.1f deg to %.1f deg)\n', bw, ...
            theta_deg(starts(main_lobe_region)), theta_deg(ends_arr(main_lobe_region)));
        bp_col = beam_pattern(:);
        if any(~above)
            fprintf('  Maximum sidelobe level: %.2f dB\n', max(bp_col(~above)));
        end
    else
        fprintf('  WARNING: could not identify main lobe region\n');
    end

    beam_linear = 10.^(beam_pattern(:) / 10);
    integral_val = trapz(theta_rad(:), beam_linear);
    DI = 10 * log10(2 * pi / integral_val);
    fprintf('  Directivity index: %.2f dB\n\n', DI);
end


% FUNCTION: BUILD_SEARCH_GRID
% Construct a consistent 2D search grid for MVDR beamforming.
% Frequency-dependent grid step (λ/N) by default; pass a positive
% grid_resolution_fixed value to override with a constant step.
function [x_scan, y_scan, X_grid, Y_grid, candidate_points, grid_res] = ...
    build_search_grid(test_freq, c_0, source_distance, array_centre_y, ...
    grid_pts_per_lambda, margin_fixed, grid_resolution_fixed)

    if isempty(grid_resolution_fixed)
        lambda = c_0 / test_freq;
        grid_resolution = lambda / grid_pts_per_lambda;
    else
        grid_resolution = grid_resolution_fixed;
    end

    max_extent = source_distance + margin_fixed;

    x_scan_points = round(2 * max_extent / grid_resolution) + 1;
    y_scan_points = round(2 * max_extent / grid_resolution) + 1;

    x_scan = linspace(-max_extent, max_extent, x_scan_points);
    y_scan = linspace(array_centre_y - max_extent, array_centre_y + max_extent, y_scan_points);

    [X_grid, Y_grid] = meshgrid(x_scan, y_scan);
    candidate_points = [X_grid(:), Y_grid(:), zeros(numel(X_grid), 1)];

    x_res = x_scan(2) - x_scan(1);
    y_res = y_scan(2) - y_scan(1);
    grid_res = mean([x_res, y_res]);
end


% FUNCTION: REFINE_PEAK_2D
% Sub-pixel peak localisation via 2D parabolic interpolation
function [est_x, est_y, refined_db] = refine_peak_2d(grid_response, x_scan, y_scan)

    % Step 1: find discrete peak
    [~, max_idx] = max(grid_response(:));
    [iy_pk, ix_pk] = ind2sub(size(grid_response), max_idx);

    est_x_discrete = x_scan(ix_pk);
    est_y_discrete = y_scan(iy_pk);

    dx = x_scan(2) - x_scan(1);
    dy = y_scan(2) - y_scan(1);

    % Step 2: check boundary — need at least one neighbour on each side
    if ix_pk < 2 || ix_pk > length(x_scan) - 1 || ...
       iy_pk < 2 || iy_pk > length(y_scan) - 1
        warning('refine_peak_2d:boundary', ...
            'Discrete peak at grid boundary — returning discrete estimate.');
        est_x = est_x_discrete;
        est_y = est_y_discrete;
        refined_db = grid_response(iy_pk, ix_pk);
        return;
    end

    % Step 3: parabolic interpolation along x at y = iy_pk
    fx_m = grid_response(iy_pk, ix_pk - 1);
    fx_0 = grid_response(iy_pk, ix_pk);
    fx_p = grid_response(iy_pk, ix_pk + 1);

    denom_x = fx_m - 2 * fx_0 + fx_p;
    if abs(denom_x) > eps
        delta_ix = (fx_m - fx_p) / (2 * denom_x);
    else
        delta_ix = 0; %flat region — no refinement possible
    end

    % Step 4: parabolic interpolation along y at x = ix_pk
    fy_m = grid_response(iy_pk - 1, ix_pk);
    fy_0 = grid_response(iy_pk, ix_pk);
    fy_p = grid_response(iy_pk + 1, ix_pk);

    denom_y = fy_m - 2 * fy_0 + fy_p;
    if abs(denom_y) > eps
        delta_iy = (fy_m - fy_p) / (2 * denom_y);
    else
        delta_iy = 0;
    end

    % Step 5: clamp offsets to ±0.5 grid cells
    delta_ix = max(-0.5, min(0.5, delta_ix));
    delta_iy = max(-0.5, min(0.5, delta_iy));

    % Step 6: convert to physical coordinates
    est_x = est_x_discrete + delta_ix * dx;
    est_y = est_y_discrete + delta_iy * dy;

    % Step 7: interpolated peak value (parabolic vertex)
    refined_db = fx_0 - (fx_m - fx_p)^2 / (8 * denom_x);

    shift_m = sqrt((delta_ix * dx)^2 + (delta_iy * dy)^2);
    fprintf('  refine_peak_2d: discrete (%.4f, %.4f) -> refined (%.4f, %.4f) | shift %.4f m\n', ...
        est_x_discrete, est_y_discrete, est_x, est_y, shift_m);
end