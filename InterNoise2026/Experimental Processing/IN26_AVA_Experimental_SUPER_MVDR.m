clc; clear all; close all;

% Two-array superimposed MVDR — experimental data (sequential recordings)
%
% Loads two independently recorded AVS array datasets, processes each
% through the AVS MVDR pipeline, then superimposes via mean-normalised
% power combination and ray intersection localisation.
%
% Two independent calibration corrections are available per array:
%   (1) GAIN CORRECTION — loaded from mic_gain_correction.mat
%       Linear gain multipliers per frequency bin, applied via interp1
%       before CSM accumulation.
%   (2) PHASE CORRECTION — defined in the calibration section below
%       Phase offsets (rad) measured at ref_freq, converted to time delays
%       and applied as exp(-1i * 2*pi*f * tau) at each bin.
%
% Written by L Marshall 29/03/2026


%% DEFINE INPUT VARIABLES %%

% Array 1 input file parameters
wav_filename_a1 = '20260328-052011(UTC)-Trial5.convent.AVA.4MIC.0.4m.4.5cmdelta.1khz-0006104060.wav';
csv_filename_a1 = 'AVA_CONVENT_Array1_Signal.csv';
convert_wav_a1 = true;

% Array 2 input file parameters
wav_filename_a2 = '20260328-063838(UTC)-Trial6.convent.AVA.4MIC.0.4m.4.5cmdelta.pos2.1khz-0010809544.wav';
csv_filename_a2 = 'AVA_CONVENT_Array2_Signal.csv';
convert_wav_a2 = true;

% Physical environment parameters
c_0 = 340; %speed of sound (m/s)
rho_0 = 1.02; %air density at STP

% Expected source position and frequency — for reference and error metrics
source_positions = [
    -0.4, 0, 0; %source 1
];
source_frequencies = [
    1000; %source 1 (Hz)
];
num_sources = size(source_positions, 1);

% Vector sensor array parameters — both arrays assumed identical configuration
N_a = 2; %number of independent arrays
N_v = 1; %vector sensors per array
d_y = 0.1; %vector sensor y-axis spacing (m)
delta = 0.04; %MEMS colocation spacing (m)

% Array 1 centre coordinates (m)
x_a1 = 0;
y_a1 = 0;
z_a1 = 0;
A1_coord = [x_a1; y_a1; z_a1];

% Array 2 centre coordinates (m) — physical location during second recording
x_a2 = 0;
y_a2 = -0.32;
z_a2 = 0;
A2_coord = [x_a2; y_a2; z_a2];

% Processing parameters
freq_limit = 2000; %maximum frequency for analysis (Hz)
overlap = 0.5; %window overlap fraction
loading = 1e-4; %regularisation parameter

% Grid parameters — frequency dependent to maintain constant points per wavelength
grid_pts_per_lambda = 20; %grid points per wavelength
margin_fixed = 0.5; %fixed physical search margin (m)

lambda = c_0 / source_frequencies(1);
grid_resolution = lambda / grid_pts_per_lambda;
max_extent = 0.4 + margin_fixed;
x_scan_points = round(2 * max_extent / grid_resolution) + 1;
y_scan_points = round(2 * max_extent / grid_resolution) + 1;


%% CALIBRATION CONFIGURATION %%

% Toggle flags — set independently per array
apply_gain_correction_a1 = true;
apply_phase_correction_a1 = false;

apply_gain_correction_a2 = true;
apply_phase_correction_a2 = false;

% Gain correction mat file — shared table assumed applicable to both arrays
gain_corr_mat_file = 'mic_gain_correction.mat';

% Phase correction — measured at ref_freq (Hz)
% Positive = channel leads, correction applies a lag
% Negative = channel lags, correction applies a lead
ref_freq = 1000;

% Array 1 phase offsets (rad at ref_freq)
phase_offset_rad_a1 = [
    +0.000000; %ch 1
    +0.216973; %ch 2
    +0.039670; %ch 3
    +0.196226; %ch 4
];

% Array 2 phase offsets (rad at ref_freq)
% Update if a separate phase measurement was taken for the array 2 position
phase_offset_rad_a2 = [
    +0.000000; %ch 1
    +0.216973; %ch 2
    +0.039670; %ch 3
    +0.196226; %ch 4
];


%% LOAD GAIN CORRECTION TABLE %%

if apply_gain_correction_a1 || apply_gain_correction_a2
    if exist(gain_corr_mat_file, 'file')
        load(gain_corr_mat_file, 'gain_corr_table', 'corr_freq_vec');
        fprintf('Loaded gain correction table: [%d freqs x %d channels]\n', ...
            size(gain_corr_table, 1), size(gain_corr_table, 2));
        fprintf('  Frequency range: %.0f - %.0f Hz\n', ...
            min(corr_freq_vec), max(corr_freq_vec));
    else
        warning('mic_gain_correction.mat not found. Gain correction disabled for both arrays.');
        apply_gain_correction_a1 = false;
        apply_gain_correction_a2 = false;
        gain_corr_table = [];
        corr_freq_vec = [];
    end
else
    gain_corr_table = [];
    corr_freq_vec = [];
end

fprintf('\n<strong>Calibration Configuration</strong>\n');
fprintf('Array 1 — Gain correction: %s | Phase correction: %s\n', ...
    string(apply_gain_correction_a1), string(apply_phase_correction_a1));
fprintf('Array 2 — Gain correction: %s | Phase correction: %s\n\n', ...
    string(apply_gain_correction_a2), string(apply_phase_correction_a2));


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

N_ma = N_v * 4; %microphones per array

% Array 1
if ~exist(csv_filename_a1, 'file')
    error('Array 1 CSV not found: %s\nRun WAV conversion first.', csv_filename_a1);
end
data_a1 = readmatrix(csv_filename_a1);
time_a1 = data_a1(:, 1);
tx_a1 = data_a1(:, 2:end).'; %rows = sensors, columns = time samples
N_a1 = size(tx_a1, 2);

if size(tx_a1, 1) ~= N_ma
    warning('Array 1: expected %d microphones but found %d in CSV.', N_ma, size(tx_a1, 1));
end
fprintf('  Array 1: %d samples from %d microphones (%.3f s)\n', ...
    N_a1, size(tx_a1, 1), time_a1(end));

% Array 2
if ~exist(csv_filename_a2, 'file')
    error('Array 2 CSV not found: %s\nRun WAV conversion first.', csv_filename_a2);
end
data_a2 = readmatrix(csv_filename_a2);
time_a2 = data_a2(:, 1);
tx_a2 = data_a2(:, 2:end).';
N_a2 = size(tx_a2, 2);

if size(tx_a2, 1) ~= N_ma
    warning('Array 2: expected %d microphones but found %d in CSV.', N_ma, size(tx_a2, 1));
end
fprintf('  Array 2: %d samples from %d microphones (%.3f s)\n\n', ...
    N_a2, size(tx_a2, 1), time_a2(end));


%% SIGNAL CONDITIONING %%

tx_a1 = tx_a1 - mean(tx_a1, 2); %remove DC offset
tx_a2 = tx_a2 - mean(tx_a2, 2);

tx_vs_arrays = {tx_a1; tx_a2};


%% ARRAY GEOMETRY SETUP %%

fprintf('<strong>Setting up array geometry</strong>\n');

array_centres = [A1_coord, A2_coord];

if N_v > 1
    array_centres(2, :) = array_centres(2, :) - (N_v - 1) * d_y / 2;
end

vs_centres_arrays = cell(N_a, 1);
mic_positions = zeros(3, N_a * N_ma);
mic_offsets = delta / 2 * [-1, 1, 1, -1; -1, -1, 1, 1; 0, 0, 0, 0];

for array = 1:N_a
    vs_centres = zeros(3, N_v);
    for vs = 1:N_v
        vs_centres(:, vs) = array_centres(:, array) + [0; (vs - 1) * d_y; 0];
        idx = (array - 1) * N_ma + (vs - 1) * 4 + (1:4);
        mic_positions(:, idx) = vs_centres(:, vs) + mic_offsets;
    end
    vs_centres_arrays{array} = vs_centres;
end

fprintf('  Array 1 centre: (%.3f, %.3f, %.3f) m\n', x_a1, y_a1, z_a1);
fprintf('  Array 2 centre: (%.3f, %.3f, %.3f) m\n', x_a2, y_a2, z_a2);
fprintf('  Array separation: %.3f m\n\n', norm(A2_coord - A1_coord));


%% FREQUENCY DOMAIN SETUP %%

% Both recordings assumed at the same sample rate — use Array 1 time vector
d_t = time_a1(2) - time_a1(1);
F_s = 1 / d_t;

fprintf('<strong>Frequency domain setup</strong>\n');
fprintf('  Sample rate: %.0f Hz | Time step: %.6f s\n', F_s, d_t);

max_binwidth = 1 / (8 * d_y * (N_ma - 1) / c_0);
size_fft = floor(F_s / max_binwidth);
fft_vec = F_s * (0:(size_fft - 1)) / size_fft;
min_freq = max_binwidth;
target_freqs = min_freq:max_binwidth:freq_limit;

bin_index = zeros(size(target_freqs));
for k = 1:length(target_freqs)
    [~, bin_index(k)] = min(abs(fft_vec - target_freqs(k)));
end
bin_index = unique(bin_index);
bin_index = bin_index(bin_index <= size_fft);
bin_freqs = fft_vec(bin_index);
num_bins = length(bin_index);

fprintf('  Using %d frequency bins (%.1f - %.1f Hz)\n\n', ...
    num_bins, min(bin_freqs), max(bin_freqs));


%% PLOT FREQUENCY SPECTRA %%

figure('Name', 'Frequency Spectra - Array 1');
f_plot = F_s * (0:floor(N_a1/2)) / N_a1;
for sensor = 1:N_ma
    Y = fft(tx_a1(sensor, :));
    Y = Y(1:floor(N_a1/2) + 1);
    Y_mag = abs(Y) / N_a1;
    Y_mag(2:end-1) = 2 * Y_mag(2:end-1);
    plot(f_plot, 20 * log10(Y_mag + eps), 'DisplayName', sprintf('Sensor %d', sensor));
    hold on;
end
hold off;
xlim([0 freq_limit]);
xlabel('Frequency (Hz)', 'FontName', 'Times New Roman', 'FontSize', 14);
ylabel('Amplitude (dB)', 'FontName', 'Times New Roman', 'FontSize', 14);
legend('show', 'Location', 'best');
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);

figure('Name', 'Frequency Spectra - Array 2');
f_plot = F_s * (0:floor(N_a2/2)) / N_a2;
for sensor = 1:N_ma
    Y = fft(tx_a2(sensor, :));
    Y = Y(1:floor(N_a2/2) + 1);
    Y_mag = abs(Y) / N_a2;
    Y_mag(2:end-1) = 2 * Y_mag(2:end-1);
    plot(f_plot, 20 * log10(Y_mag + eps), 'DisplayName', sprintf('Sensor %d', sensor));
    hold on;
end
hold off;
xlim([0 freq_limit]);
xlabel('Frequency (Hz)', 'FontName', 'Times New Roman', 'FontSize', 14);
ylabel('Amplitude (dB)', 'FontName', 'Times New Roman', 'FontSize', 14);
legend('show', 'Location', 'best');
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);


%% CREATE AVS SNAPSHOTS AND CROSS-SPECTRAL MATRICES %%

fprintf('<strong>Creating snapshots and cross-spectral matrices</strong>\n');

window = hanning(size_fft)';

% Package calibration settings per array for passing into create_vs_csm
apply_gain = {apply_gain_correction_a1, apply_gain_correction_a2};
apply_phase = {apply_phase_correction_a1, apply_phase_correction_a2};
phase_offsets = {phase_offset_rad_a1, phase_offset_rad_a2};

snapshots_vs_arrays = cell(N_a, 1);
r_vs_arrays = cell(N_a, 1);

for array = 1:N_a
    snapshots_vs_arrays{array} = make_snapshots(tx_vs_arrays{array}, size_fft, overlap, window);
    fprintf('  Array %d: created %d snapshots\n', array, size(snapshots_vs_arrays{array}, 3));

    r_vs_arrays{array} = create_vs_csm(snapshots_vs_arrays{array}, bin_index, delta, ...
        bin_freqs, rho_0, N_v, num_bins, c_0, ...
        apply_gain{array}, gain_corr_table, corr_freq_vec, ...
        apply_phase{array}, phase_offsets{array}, ref_freq);

    fprintf('  Array %d: created %dx%dx%d cross-spectral matrix\n', ...
        array, size(r_vs_arrays{array}));
end
fprintf('\n');


%% DEFINE GRID SEARCH AREA %%

fprintf('<strong>Defining search area</strong>\n');

system_centre_y = -0.16; %geometric midpoint of full aperture (m)
[x_scan, y_scan, X_grid, Y_grid, candidate_points, grid_res] = ...
    build_search_grid(source_frequencies(1), c_0, 0.4, ...
    system_centre_y, grid_pts_per_lambda, margin_fixed);

fprintf('  Search grid: %d x %d points, resolution: %.4f m\n\n', ...
    length(x_scan), length(y_scan), grid_res);


%% PERFORM MVDR BEAMFORMING %%

fprintf('<strong>Running MVDR Beamforming</strong>\n');
fprintf('  Regularisation parameter: %.2e\n\n', loading);

response_db_arrays = cell(N_a, 1);

for array = 1:N_a
    fprintf('  Processing Array %d...\n', array);
    response_db_arrays{array} = mvdr_vs_beamforming(r_vs_arrays{array}, ...
        vs_centres_arrays{array}, candidate_points, bin_freqs, c_0, rho_0, ...
        loading, num_bins, N_v);
    fprintf('    Array %d response range: [%.2f, %.2f] dB\n', array, ...
        min(response_db_arrays{array}), max(response_db_arrays{array}));
end


%% PLOT INDIVIDUAL ARRAY 2D SCANS %%

fprintf('\n2D Scans - Individual Arrays:\n');
for array = 1:N_a
    fprintf('  Plotting Array %d response...\n', array);
    plot_2dscan_mvdr(r_vs_arrays{array}, vs_centres_arrays{array}, ...
        candidate_points, bin_freqs, c_0, rho_0, y_scan, x_scan, ...
        X_grid, Y_grid, loading, num_bins, N_v, num_sources, source_positions, ...
        response_db_arrays{array}, array, grid_res, vs_centres_arrays, N_a);
end


%% COMBINE AND SUPERIMPOSE ARRAY RESPONSES %%

fprintf('\n<strong>Combining array responses</strong>\n');

% Product of linear responses before normalisation
response_linear_combined = 1;
for k = 1:N_a
    response_linear_combined = response_linear_combined .* 10.^(response_db_arrays{k} / 10);
end
response_db_combined = 10 * log10(response_linear_combined + eps);
response_db_combined = response_db_combined - max(response_db_combined);

fprintf('  Combined response range: [%.2f, %.2f] dB\n', ...
    min(response_db_combined), max(response_db_combined));
fprintf('  Plotting superimposed response and localising sources...\n');

[est_positions_combined, ~, response_db_norm_combined] = ...
    superimposed_localisation(response_db_combined, X_grid, Y_grid, y_scan, x_scan, ...
    num_sources, source_positions, N_a, grid_res, response_db_arrays, vs_centres_arrays);


%% PERFORMANCE METRICS %%

fprintf('\n<strong>PERFORMANCE METRICS</strong>\n');

response_names = [arrayfun(@(i) sprintf('Array %d', i), 1:N_a, 'UniformOutput', false), ...
    'Superimposed'];
response_db_arrays_norm = cellfun(@(r) r - max(r), response_db_arrays, 'UniformOutput', false);
responses_to_analyse = [response_db_arrays_norm', {response_db_norm_combined}];

% Angular reference: individual array centres for per-array cases,
% system midpoint for the superimposed case
array_geometric_centres = zeros(2, N_a);
for array = 1:N_a
    array_geometric_centres(:, array) = mean(vs_centres_arrays{array}(1:2,:), 2);
end
angular_refs_individual = mat2cell(array_geometric_centres, 2, ones(1, N_a));
angular_refs = [angular_refs_individual, {mean(array_geometric_centres, 2)}];

num_cases = length(responses_to_analyse);
all_cases_results = cell(num_cases, 1);
min_separation = max([0.25, 3 * grid_res]);

for case_idx = 1:num_cases
    fprintf('\n<strong>%s:</strong>\n', response_names{case_idx});

    grid_response = reshape(responses_to_analyse{case_idx}, length(y_scan), length(x_scan));
    ref_pos = angular_refs{case_idx};

    results = struct('est_positions', zeros(num_sources, 2), ...
        'true_positions', zeros(num_sources, 2), ...
        'errors', zeros(num_sources, 2), ...
        'MSE', zeros(num_sources, 1), ...
        'MSE_percent', zeros(num_sources, 1), ...
        'radial_error', zeros(num_sources, 1), ...
        'angular_error', zeros(num_sources, 1), ...
        'peak_response', zeros(num_sources, 1));

    est_positions = zeros(num_sources, 2);
    peak_responses = zeros(num_sources, 1);

    % Superimposed case: use ray intersection estimate directly
    % Individual array cases: standard peak search on normalised response
    if case_idx == num_cases
        est_positions = est_positions_combined;
        for src = 1:num_sources
            y_idx = find(abs(Y_grid(:,1) - est_positions(src,2)) == ...
                min(abs(Y_grid(:,1) - est_positions(src,2))), 1);
            x_idx = find(abs(X_grid(1,:) - est_positions(src,1)) == ...
                min(abs(X_grid(1,:) - est_positions(src,1))), 1);
            peak_responses(src) = grid_response(y_idx, x_idx);
        end
    else
        response_work = grid_response;
        for src = 1:num_sources
            [peak_x, peak_y, peak_responses(src)] = refine_peak_2d(response_work, x_scan, y_scan);
            est_positions(src, :) = [peak_x, peak_y];
            if num_sources > 1 && src < num_sources
                distances = sqrt((X_grid - est_positions(src,1)).^2 + ...
                    (Y_grid - est_positions(src,2)).^2);
                response_work(distances < min_separation) = -Inf;
            end
        end
    end

    for src = 1:num_sources
        fprintf('  <strong>SOURCE %d: (%.3f, %.3f) m @ %.0f Hz</strong>\n', ...
            src, source_positions(src, 1:2), source_frequencies(src));

        true_pos = source_positions(src, 1:2);
        est_pos = est_positions(src, :);
        errors = true_pos - est_pos;
        radial_error = norm(errors);
        MSE = sum(errors.^2);
        MSE_percent = 100 * MSE / (grid_res^2);

        true_angle = atan2(true_pos(2) - ref_pos(2), true_pos(1) - ref_pos(1));
        est_angle = atan2(est_pos(2) - ref_pos(2), est_pos(1) - ref_pos(1));
        angular_error_deg = rad2deg(abs(true_angle - est_angle));
        if angular_error_deg > 180
            angular_error_deg = 360 - angular_error_deg;
        end

        results.est_positions(src, :) = est_pos;
        results.true_positions(src, :) = true_pos;
        results.errors(src, :) = errors;
        results.MSE(src) = MSE;
        results.MSE_percent(src) = MSE_percent;
        results.radial_error(src) = radial_error;
        results.angular_error(src) = angular_error_deg;
        results.peak_response(src) = peak_responses(src);

        fprintf('  Estimated position: (%.4f, %.4f) m\n', est_pos);
        fprintf('  Position error: dx = %.4f m, dy = %.4f m\n', errors);
        fprintf('  Radial error: %.4f m\n', radial_error);
        fprintf('  Mean squared error: %.4f m^2\n', MSE);
        fprintf('  MSE (normalised): %.3f%% of grid resolution\n', MSE_percent);
        fprintf('  Angular error: %.3f deg\n', angular_error_deg);
        fprintf('  Peak response: %.3f dB\n\n', peak_responses(src));
    end

    all_cases_results{case_idx} = results;
end

% Comparative summary
fprintf('\n<strong>COMPARATIVE ANALYSIS:</strong>\n\n');

fprintf('Radial Error Comparison:\n');
fprintf('%-10s ', 'Source');
for c = 1:num_cases; fprintf('%-15s ', response_names{c}); end
fprintf('\n%s\n', repmat('-', 1, 10 + 15 * num_cases));
for src = 1:num_sources
    fprintf('%-10d ', src);
    for c = 1:num_cases
        fprintf('%-15.4f ', all_cases_results{c}.radial_error(src));
    end
    fprintf('\n');
end

fprintf('\nAngular Error Comparison:\n');
fprintf('%-10s ', 'Source');
for c = 1:num_cases; fprintf('%-15s ', response_names{c}); end
fprintf('\n%s\n', repmat('-', 1, 10 + 15 * num_cases));
for src = 1:num_sources
    fprintf('%-10d ', src);
    for c = 1:num_cases
        fprintf('%-15.3f ', all_cases_results{c}.angular_error(src));
    end
    fprintf('\n');
end
fprintf('\n');

all_results = all_cases_results{end};


%% BEAM PATTERN ANALYSIS %%

fprintf('\n<strong>GENERATING BEAM PATTERNS</strong>\n');

distances_all = [];
for array = 1:N_a
    array_centre_2d = mean(vs_centres_arrays{array}(1:2,:), 2);
    distances = sqrt(sum((source_positions(:,1:2) - array_centre_2d').^2, 2));
    distances_all = [distances_all; distances];
end
radius = median(distances_all);
fprintf('Using radius = %.3f m (median source distance across all arrays)\n', radius);

for array = 1:N_a
    array_centre_2d = mean(vs_centres_arrays{array}(1:2,:), 2);
    [~, ~, ~] = compute_nearfield_beam_pattern(...
        r_vs_arrays{array}, vs_centres_arrays{array}, array_centre_2d, ...
        source_positions(:,1:2), bin_freqs, c_0, rho_0, loading, num_bins, N_v, ...
        radius, sprintf('Array %d Near-Field Beam Pattern', array));
end


%% FUNCTION DEFINITIONS %%


% FUNCTION: MAKE_SNAPSHOTS
% Make time-domain snapshots of input signal
function snapshots = make_snapshots(tx_vs, size_fft, overlap, window)

    step = round(overlap * size_fft);
    num_snap = floor((size(tx_vs, 2) - size_fft) / step) + 1;

    start_idx = 1 + (0:(num_snap - 1)) * step;
    idx_matrix = start_idx + (0:size_fft - 1)';

    snapshots = tx_vs(:, idx_matrix);
    snapshots = reshape(snapshots, size(tx_vs, 1), size_fft, num_snap);
    snapshots = snapshots .* reshape(window, [1, size_fft, 1]);
end


% FUNCTION: CREATE_VS_CSM
% Convert 4-mic pressure snapshots to AVS [P, Vx, Vy] outputs and accumulate
% the cross-spectral matrix per frequency bin, with optional gain and phase
% correction applied in the frequency domain before CSM accumulation.
%
% Gain correction: gain_corr_table interpolated at each bin frequency and
% applied column-wise to frequency-domain snapshots.
%
% Phase correction: phase_offset_rad (measured at ref_freq) converted to
% time delay tau = phi / (2*pi*ref_freq), applied as exp(-1i*2*pi*f*tau)
% at each bin — correctly scales phase angle with frequency.
function r_vs = create_vs_csm(snapshots_vs, bin_index, delta, bin_freqs, rho_0, N_v, num_bins, c_0, ...
    apply_gain, gain_corr_table, corr_freq_vec, apply_phase, phase_offset_rad, ref_freq)

    signal_dim = 3 * N_v;
    fx = fft(snapshots_vs, [], 2);
    fx = fx(:, bin_index, :);
    r_vs = zeros(signal_dim, signal_dim, num_bins);
    num_snap = size(snapshots_vs, 3);

    for jf = 1:num_bins
        freq = bin_freqs(jf);
        fx_bin = squeeze(fx(:, jf, :)); %(N_ch x num_snap)

        % Gain correction — interpolated from table at current bin frequency
        if apply_gain && ~isempty(gain_corr_table)
            gain_vec = interp1(corr_freq_vec, gain_corr_table, freq, 'linear', 'extrap');
            fx_bin = fx_bin .* gain_vec(:);
        end

        % Phase correction — time-delay model scaled to bin frequency
        if apply_phase && ~isempty(phase_offset_rad)
            tau_vec = phase_offset_rad(:) / (2 * pi * ref_freq);
            phase_corr_vec = exp(-1i * 2 * pi * freq * tau_vec);
            fx_bin = fx_bin .* phase_corr_vec;
        end

        R_freq = zeros(signal_dim, signal_dim);

        for snap = 1:num_snap
            fx_snap = fx_bin(:, snap);
            [p_vs, vx_vs, vy_vs] = vs_outputs(fx_snap, delta, freq, rho_0, N_v, c_0);

            rho_c = rho_0 * c_0;
            vx_vs = rho_c * vx_vs;
            vy_vs = rho_c * vy_vs;

            vs_signal = zeros(signal_dim, 1);
            for n = 1:N_v
                idx = (n - 1) * 3 + (1:3);
                vs_signal(idx) = [p_vs(n); vx_vs(n); vy_vs(n)];
            end

            R_freq = R_freq + (vs_signal * vs_signal');
        end

        r_vs(:,:,jf) = R_freq / num_snap;
    end
end


% FUNCTION: MVDR_VS_BEAMFORMING
% MVDR beamforming for vector sensor array with adaptive diagonal loading
function response_db = mvdr_vs_beamforming(r_vs, vs_centres, candidate_points, ...
    bin_freqs, c_0, rho_0, loading, num_bins, N_v)

    num_points = size(candidate_points, 1);
    mvdr_responses = zeros(num_points, num_bins);

    for jf = 1:num_bins
        freq = bin_freqs(jf);
        rf = squeeze(r_vs(:,:,jf));

        [ur, er] = eig(rf);
        erv = real(diag(er));
        erv = max(erv, 1e-12);
        lambda = loading * max(erv);

        for n = 1:num_points
            source_pos = candidate_points(n,:).';
            v_vs = vs_steering_vector(vs_centres, source_pos, freq, c_0, rho_0, N_v);

            rxv = (ur * diag(1 ./ (erv + lambda)) * ur') * v_vs;
            denominator = v_vs' * rxv;

            if abs(denominator) > 1e-12
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
% Returns estimated source positions and peak response value
function [all_est_positions, scan_2d_max_db] = plot_2dscan_mvdr(...
    r_vs, vs_centres, candidate_points, bin_freqs, c_0, rho_0, y_scan, x_scan, ...
    X_grid, Y_grid, loading, num_bins, N_v, num_sources, source_positions, response_db, ...
    array_num, grid_res, vs_centres_arrays, N_a)

    if nargin < 16
        response_db = mvdr_vs_beamforming(r_vs, vs_centres, candidate_points, ...
            bin_freqs, c_0, rho_0, loading, num_bins, N_v);
        array_num = 1;
    elseif nargin < 17
        array_num = 1;
    end

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

    [scan_2d_max_db, ~] = max(grid_response(:));
    min_separation = max([0.25, 3 * grid_res]);
    all_est_positions = zeros(num_sources, 2);
    response_work = grid_response;

    for src = 1:num_sources
        [~, max_idx] = max(response_work(:));
        [y_idx, x_idx] = ind2sub(size(grid_response), max_idx);
        all_est_positions(src, :) = [X_grid(y_idx, x_idx), Y_grid(y_idx, x_idx)];
        if num_sources > 1 && src < num_sources
            distances = sqrt((X_grid - all_est_positions(src,1)).^2 + ...
                (Y_grid - all_est_positions(src,2)).^2);
            response_work(distances < min_separation) = -Inf;
        end
    end

    % True source — black cross
    plot(source_positions(:,1), source_positions(:,2), 'kx', ...
        'MarkerSize', 5, 'LineWidth', 1.5, 'DisplayName', 'True source');

    % Estimated source — black hollow circle
    plot(all_est_positions(:,1), all_est_positions(:,2), 'ko', ...
        'MarkerSize', 10, 'LineWidth', 1.5, 'MarkerFaceColor', 'none', ...
        'DisplayName', 'Estimated source');

    % Array centre — filled white square for the array being plotted
    ac = mean(vs_centres_arrays{array_num}(1:2,:), 2);
    plot(ac(1), ac(2), 'ks', 'MarkerSize', 9, 'LineWidth', 1, ...
        'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'DisplayName', 'Array');

    legend('Location', 'northeast', 'FontName', 'Times New Roman', ...
        'FontSize', 10, 'Color', 'w', 'TextColor', 'k', 'Box', 'on');
end


% FUNCTION: SUPERIMPOSED_LOCALISATION
% Mean-normalised power combination with ray intersection as source estimate
%
% Stage 1: each array response normalised to unit mean power before
% multiplication — reflects regions of high response in both arrays
% regardless of absolute power differences between recordings.
%
% Stage 2: geometric intersection of the two bearing rays taken as
% the estimated source position.
function [all_est_positions, max_db, response_db_norm_combined_out] = ...
    superimposed_localisation(response_db_combined, X_grid, Y_grid, y_scan, x_scan, ...
    num_sources, source_positions, N_a, grid_res, response_db_arrays, vs_centres_arrays)

    min_separation = max([0.25, 3 * grid_res]);

    % Stage 1: mean-normalised power combination
    if nargin >= 11 && ~isempty(response_db_arrays)
        response_linear_combined_norm = ones(size(response_db_arrays{1}));
        for k = 1:N_a
            linear_k = 10.^(response_db_arrays{k} / 10);
            linear_k_norm = linear_k / mean(linear_k);
            response_linear_combined_norm = response_linear_combined_norm .* linear_k_norm;
        end
        response_db_norm_combined = 10 * log10(response_linear_combined_norm + eps);
        response_db_norm_combined = response_db_norm_combined - max(response_db_norm_combined);
    else
        response_db_norm_combined = response_db_combined;
    end

    grid_response_norm = reshape(response_db_norm_combined, length(y_scan), length(x_scan));
    grid_response_display = reshape(response_db_combined, length(y_scan), length(x_scan));
    [max_db, ~] = max(grid_response_display(:));

    response_db_arrays_norm = cellfun(@(r) r - max(r), response_db_arrays, 'UniformOutput', false);

    % Stage 2: ray intersection estimate from per-array peak bearings
    % Solve p1 + t1*d1 = p2 + t2*d2 for the 2D intersection point
    ray_intersection = [];

    if N_a >= 2 && ~isempty(vs_centres_arrays)
        ray_origins = zeros(2, N_a);
        ray_directions = zeros(2, N_a);
        ray_valid = false(1, N_a);

        for array = 1:N_a
            array_centre_2d = mean(vs_centres_arrays{array}(1:2,:), 2);
            ray_origins(:, array) = array_centre_2d;

            arr_response = reshape(response_db_arrays_norm{array}, length(y_scan), length(x_scan));
            [peak_x, peak_y, ~] = refine_peak_2d(arr_response, x_scan, y_scan);
            peak_pos = [peak_x; peak_y];

            dir = peak_pos - array_centre_2d;
            dir_norm = norm(dir);
            if dir_norm > 1e-6
                ray_directions(:, array) = dir / dir_norm;
                ray_valid(array) = true;
                fprintf('  Array %d: peak at (%.3f, %.3f) m | bearing = %.1f deg\n', ...
                    array, peak_pos(1), peak_pos(2), rad2deg(atan2(dir(1), dir(2))));
            end
        end

        valid_idx = find(ray_valid);
        if length(valid_idx) >= 2
            i1 = valid_idx(1); i2 = valid_idx(2);
            p1 = ray_origins(:,i1); d1 = ray_directions(:,i1);
            p2 = ray_origins(:,i2); d2 = ray_directions(:,i2);
            A_mat = [d1, -d2];
            b_vec = p2 - p1;

            if abs(det(A_mat)) > 1e-6
                t_params = A_mat \ b_vec;
                intersection_2d = p1 + t_params(1) * d1;

                if intersection_2d(1) >= min(x_scan) && intersection_2d(1) <= max(x_scan) && ...
                        intersection_2d(2) >= min(y_scan) && intersection_2d(2) <= max(y_scan)
                    ray_intersection = intersection_2d;
                    fprintf('  Ray intersection: (%.3f, %.3f) m\n', ...
                        ray_intersection(1), ray_intersection(2));
                else
                    fprintf('  Ray intersection (%.3f, %.3f) m outside grid — falling back to global peak\n', ...
                        intersection_2d(1), intersection_2d(2));
                end
            else
                fprintf('  WARNING: array rays are parallel — falling back to global peak\n');
            end
        end
    end

    % Source position estimate
    all_est_positions = zeros(num_sources, 2);
    for src = 1:num_sources
        if src == 1 && ~isempty(ray_intersection)
            all_est_positions(src, :) = ray_intersection';
            fprintf('  Source %d: estimated from ray intersection\n', src);
        else
            [~, max_idx] = max(grid_response_norm(:));
            [y_idx, x_idx] = ind2sub(size(grid_response_norm), max_idx);
            all_est_positions(src, :) = [X_grid(y_idx, x_idx), Y_grid(y_idx, x_idx)];
            fprintf('  Source %d: no valid intersection — estimated from global peak\n', src);
        end
        fprintf('    Source %d: estimated at (%.4f, %.4f) m\n', src, ...
            all_est_positions(src,1), all_est_positions(src,2));
    end

    response_db_norm_combined_out = response_db_norm_combined;

    % Publication style plot
    figure('Color', 'w', 'Position', [100 100 520 440]);
    imagesc(x_scan, y_scan, grid_response_norm);
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

    % True source — black cross
    plot(source_positions(:,1), source_positions(:,2), 'kx', ...
        'MarkerSize', 5, 'LineWidth', 1.5, 'DisplayName', 'True source');

    % Estimated source — black hollow circle
    plot(all_est_positions(:,1), all_est_positions(:,2), 'ko', ...
        'MarkerSize', 10, 'LineWidth', 1.5, 'MarkerFaceColor', 'none', ...
        'DisplayName', 'Estimated source');

    % Array centres — filled white squares, first only labelled
    for array = 1:N_a
        ac = mean(vs_centres_arrays{array}(1:2,:), 2);
        if array == 1
            plot(ac(1), ac(2), 'ks', 'MarkerSize', 9, 'LineWidth', 1, ...
                'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'DisplayName', 'Array');
        else
            plot(ac(1), ac(2), 'ks', 'MarkerSize', 9, 'LineWidth', 1, ...
                'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');
        end
    end

    legend('Location', 'northeast', 'FontName', 'Times New Roman', ...
        'FontSize', 10, 'Color', 'w', 'TextColor', 'k', 'Box', 'on');
end


% FUNCTION: VS_OUTPUTS
% Vector sensor pressure and velocity outputs via least-squares gradient estimation
function [p_vs, vx_vs, vy_vs] = vs_outputs(tx_freq, delta, freq, rho_0, N_v, c_0)

    omega = 2 * pi * freq;
    p_vs = zeros(N_v, 1);
    vx_vs = zeros(N_v, 1);
    vy_vs = zeros(N_v, 1);

    for vs = 1:N_v
        idx = (vs - 1) * 4 + (1:4);

        mic_pos = [-delta/2, -delta/2;
                    delta/2, -delta/2;
                    delta/2,  delta/2;
                   -delta/2,  delta/2];

        M = [ones(4, 1), mic_pos];
        p_vector = tx_freq(idx);
        coeffs = M \ p_vector;

        p_vs(vs) = coeffs(1);
        vx_vs(vs) = -coeffs(2) / (1i * omega * rho_0);
        vy_vs(vs) = -coeffs(3) / (1i * omega * rho_0);
    end
end


% FUNCTION: VS_STEERING_VECTOR
% Near-field steering vector for acoustic vector sensor array
function v_vs = vs_steering_vector(vs_centres, source_pos, freq, c_0, rho_0, N_v)

    k = 2 * pi * freq / c_0;
    omega = 2 * pi * freq;
    v_vs = zeros(3 * N_v, 1);

    for n = 1:N_v
        r_vec = vs_centres(:,n) - source_pos;
        r = norm(r_vec);
        r_hat = r_vec / r;

        p_steer = exp(-1i * k * r) / r;
        common_term = exp(-1i * k * r) / (1i * omega * rho_0 * r^2);
        rho_c = rho_0 * c_0;
        vx_steer = rho_c * common_term * (1 + 1i * k * r) * r_hat(1);
        vy_steer = rho_c * common_term * (1 + 1i * k * r) * r_hat(2);

        idx = (n - 1) * 3 + (1:3);
        v_vs(idx) = [p_steer; vx_steer; vy_steer];
    end
end


% FUNCTION: COMPUTE_NEARFIELD_BEAM_PATTERN
% Evaluate MVDR response at each azimuth at fixed range radius
% Returns beam pattern statistics to command window
function [theta_deg, beam_pattern, fig_handle] = compute_nearfield_beam_pattern(...
    r_vs, vs_centres, array_centre, source_positions, bin_freqs, c_0, rho_0, ...
    loading, num_bins, N_v, radius, plot_title)

    fprintf('\n<strong>Computing Near-Field Beam Pattern:</strong>\n');
    fprintf('  Range: %.3f m\n', radius);
    fprintf('  Array centre: (%.3f, %.3f) m\n', array_centre);

    num_angles = 360;
    theta_deg = linspace(0, 360, num_angles + 1);
    theta_deg = theta_deg(1:end-1); %remove duplicate at 360 deg
    theta_rad = deg2rad(theta_deg);
    mvdr_responses = zeros(num_angles, num_bins);

    for jf = 1:num_bins
        freq = bin_freqs(jf);
        rf = squeeze(r_vs(:,:,jf));

        [ur, er] = eig(rf);
        erv = real(diag(er));
        erv = max(erv, 1e-12);
        lambda = loading * max(erv);

        for angle_idx = 1:num_angles
            % 0 deg is along +Y axis (north), increases clockwise
            x_pos = array_centre(1) + radius * sin(theta_rad(angle_idx));
            y_pos = array_centre(2) + radius * cos(theta_rad(angle_idx));
            source_pos = [x_pos; y_pos; 0];

            v_vs = vs_steering_vector(vs_centres, source_pos, freq, c_0, rho_0, N_v);
            rxv = (ur * diag(1 ./ (erv + lambda)) * ur') * v_vs;
            denominator = v_vs' * rxv;

            if abs(denominator) > 1e-12
                vmvdr = rxv / denominator;
                mvdr_responses(angle_idx, jf) = abs(vmvdr' * rf * vmvdr);
            else
                mvdr_responses(angle_idx, jf) = 0;
            end
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
    ax.ThetaZeroLocation = 'top';
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
% Construct a consistent 2D search grid for MVDR beamforming
% Written by L Marshall — grid consistency fix for InterNoise 2026 scripts
function [x_scan, y_scan, X_grid, Y_grid, candidate_points, grid_res] = ...
    build_search_grid(test_freq, c_0, source_distance, array_centre_y, ...
    grid_pts_per_lambda, margin_fixed)

    lambda = c_0 / test_freq;
    grid_resolution = lambda / grid_pts_per_lambda;

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
% Written by L Marshall — grid quantisation fix for InterNoise 2026 scripts
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