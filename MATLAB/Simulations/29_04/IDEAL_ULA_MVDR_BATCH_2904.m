clc; clear all; close all;

% Batch frequency sweep testing — ULA MVDR
% Runs signal generation and ULA MVDR processing at a fixed source distance
% across a range of test frequencies for localisation performance comparison
%
% Sweep frequencies: 630, 800, 1000, 1250, 1600, 2000, 2500 Hz
% Fixed source distance: 0.4 m at 180 deg from +X axis
%
% For each frequency this script:
%   (1) Generates a synthetic signal via ULA_OUTPUT_BATCH.m
%   (2) Runs MVDR beamforming and extracts the estimated source position
%   (3) Extracts radial error, angular error, and beam pattern metrics
%   (4) Computes MVDR response along the ray from array centre through the
%       estimated source position for radial localisation assessment
%   (5) Saves all results to CSV and .mat for post-processing
%
% Written by L Marshall 29/03/2026
% Updated for ideal 8-mic ULA by L Marshall 17/04/2026
% 1D ray response analysis added by L Marshall 17/04/2026
% Corrected FFT bin selection by L Marshall 29/04/2026
% Fixed mic ordering convention by L Marshall 29/04/2026
%
% Changes (29/04/2026):
%   - CSV reader updated for single-signal-per-channel format (was reading
%     monopole + sinusoidal columns and summing them — those two components
%     were redundant at one frequency)
%   - FFT length and bin selection now centred on the test frequency rather
%     than driven by the spatial-resolution lower-bound formula. Analysis
%     bandwidth is a narrow window around f_test where the signal energy
%     actually lives — appropriate for performance characterisation of a
%     known monochromatic test source. The previous broadband scheme
%     integrated many noise-only bins, biasing low-frequency results
%   - Mic ordering in MVDR script now matches the OUTPUT script (descending
%     y, M1 at the top). Previously the MVDR script used ascending y while
%     the OUTPUT script wrote in descending y, producing a steering vector
%     mirror-reflected through the array centre. This was masked by the
%     broadband bin scheme (which smeared the peak across the bearing
%     direction so the mirror image blended with the true peak — visible
%     as the broad lobe in paper Figure 7a) but became apparent once the
%     analysis was concentrated on the signal frequency


%% BATCH TEST PARAMETERS %%

% Test frequencies (Hz)
batch_test_frequencies = [1000];
batch_num_tests = length(batch_test_frequencies);

% Fixed source geometry
source_distance_m = 0.4; %fixed source distance (m)
source_angle_deg = 180; %source angle from +X axis (deg)
c_0 = 340; %speed of sound (m/s)

% Array parameters — must match ULA_OUTPUT_BATCH.m
N_mics = 8; %microphones in array
x_a = 0; %array centre x coordinate (m)
y_a = -0.16; %array centre y coordinate (m)
z_a = 0; %array centre z coordinate (m)
total_aperture = 0.36; %total array aperture (m)
d_y = total_aperture / 7; %inter-element spacing (m) — 0.051429 m

% Processing parameters
overlap = 0.5; %window overlap fraction
loading = 1e-4; %regularisation parameter
analysis_halfwidth_hz = 10; %half-width of analysis band around f_test (Hz)
target_binwidth_hz = 1; %target FFT bin spacing within analysis band (Hz)

% Grid parameters — frequency-dependent by default to maintain constant
% spatial resolution relative to wavelength, with a fixed physical search
% margin. Set grid_resolution_fixed to a positive value (e.g. 0.005) to
% override with a constant grid step across all frequencies — useful for
% isolating grid quantisation effects from genuine frequency-dependent
% performance variation.
grid_pts_per_lambda = 20; %grid points per wavelength at each test frequency
grid_resolution_fixed = []; %override grid step (m); [] = use λ/grid_pts_per_lambda
% grid_resolution_fixed = 0.005;
margin_fixed = 0.5; %fixed physical margin around source distance (m)

% 1D ray analysis parameters
ray_num_points = 200; %number of radial sample points along ray
ray_range_lambdas = [0.1, 5.0]; %radial extent of ray in units of lambda

% Create results folder
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
results_folder = sprintf('freq_sweep_ULA_%s', timestamp);
mkdir(results_folder);

fprintf('\n<strong>BATCH FREQUENCY SWEEP TEST - ULA MVDR</strong>\n');
fprintf('Testing %d frequencies: ', batch_num_tests);
fprintf('%.0f  ', batch_test_frequencies);
fprintf('Hz\n');
fprintf('Fixed source distance: %.3f m at %.0f deg\n', source_distance_m, source_angle_deg);
fprintf('N_mics: %d | d_y: %.5f m\n', N_mics, d_y);
if isempty(grid_resolution_fixed)
    fprintf('Grid: %d pts/lambda (frequency-dependent) | Margin: %.2f m\n', ...
        grid_pts_per_lambda, margin_fixed);
else
    fprintf('Grid: %.4f m fixed | Margin: %.2f m\n', ...
        grid_resolution_fixed, margin_fixed);
end
fprintf('Ray: %d points | Range: [%.2f, %.2f] lambda\n', ...
    ray_num_points, ray_range_lambdas(1), ray_range_lambdas(2));
fprintf('Results folder: %s\n\n', results_folder);


%% INITIALISE RESULTS STORAGE %%

batch_res_frequency_hz = zeros(batch_num_tests, 1);
batch_res_lambda_m = zeros(batch_num_tests, 1);
batch_res_source_x = zeros(batch_num_tests, 1);
batch_res_source_y = zeros(batch_num_tests, 1);
batch_res_grid_resolution = zeros(batch_num_tests, 1);
batch_res_radial_error = zeros(batch_num_tests, 1);
batch_res_angular_error = zeros(batch_num_tests, 1);
batch_res_beamwidth = zeros(batch_num_tests, 1);
batch_res_sidelobe = zeros(batch_num_tests, 1);
batch_res_DI = zeros(batch_num_tests, 1);
batch_res_ray_peak = zeros(batch_num_tests, 1);


%% RUN BATCH TESTS %%

for batch_test_idx = 1:batch_num_tests

    fprintf('\n========================================\n');
    fprintf('<strong>TEST %d/%d: f = %.0f Hz</strong>\n', ...
        batch_test_idx, batch_num_tests, batch_test_frequencies(batch_test_idx));
    fprintf('========================================\n');

    batch_test_freq = batch_test_frequencies(batch_test_idx);
    lambda = c_0 / batch_test_freq;

    % Source position at fixed distance and angle
    source_x = source_distance_m * cosd(source_angle_deg);
    source_y = source_distance_m * sind(source_angle_deg);

    % Grid step — frequency-dependent (λ/N) by default, or the override
    % value if grid_resolution_fixed is set. The actual search grid is
    % constructed inside build_search_grid() using this same value
    if isempty(grid_resolution_fixed)
        grid_resolution = lambda / grid_pts_per_lambda;
    else
        grid_resolution = grid_resolution_fixed;
    end

    % Analysis band — narrow window centred on the test frequency,
    % matching the known monochromatic source spectrum
    f_band_lo = batch_test_freq - analysis_halfwidth_hz;
    f_band_hi = batch_test_freq + analysis_halfwidth_hz;

    % Store frequency-dependent parameters
    batch_res_frequency_hz(batch_test_idx) = batch_test_freq;
    batch_res_lambda_m(batch_test_idx) = lambda;
    batch_res_source_x(batch_test_idx) = source_x;
    batch_res_source_y(batch_test_idx) = source_y;
    batch_res_grid_resolution(batch_test_idx) = grid_resolution;

    fprintf('Source position: (%.3f, %.3f) m\n', source_x, source_y);
    fprintf('Lambda: %.4f m | Grid resolution: %.4f m\n', lambda, grid_resolution);
    fprintf('Analysis band: %.0f-%.0f Hz\n', f_band_lo, f_band_hi);


    % STEP 1: CLEAR PREVIOUS RUN VARIABLES %

    vars_to_clear = {'source_positions', 'source_frequencies'};
    for v = 1:length(vars_to_clear)
        if exist(vars_to_clear{v}, 'var')
            clear(vars_to_clear{v});
        end
    end


    % STEP 2: RUN SIGNAL GENERATION SCRIPT %

    fprintf('\nGenerating signal at r=%.3f m, f=%.0f Hz...\n', ...
        source_distance_m, batch_test_freq);

    batch_test_source_x = source_x;
    batch_test_source_y = source_y;
    batch_csv_name = fullfile(results_folder, ...
        sprintf('signal_%02d_%.0fhz.csv', batch_test_idx, batch_test_freq));

    fprintf('  batch_test_freq = %.0f Hz\n', batch_test_freq);
    fprintf('  batch_test_source_x = %.3f m\n', batch_test_source_x);
    fprintf('  batch_test_source_y = %.3f m\n', batch_test_source_y);
    fprintf('  batch_csv_name = %s\n', batch_csv_name);

    run('IDEAL_ULA_OUTPUT_BATCH_2904.m');
    fprintf('  Signal generation complete\n');


    % STEP 3: LOAD SIGNAL FROM CSV %
    % CSV format: [time, mic_1, mic_2, ..., mic_Nr] — single signal per channel

    fprintf('\nLoading signal from: %s\n', csv_filename);

    if ~exist(csv_filename, 'file')
        error('Signal file not found: %s', csv_filename);
    end

    data = readmatrix(csv_filename);
    time = data(:, 1);
    [N, numCols] = size(data);
    Nr = numCols - 1; %number of receiving elements

    fprintf('  Loaded %d samples from %d microphones\n', N, Nr);


    % STEP 4: ASSEMBLE SIGNAL MATRIX %
    % Rows = sensors, columns = time samples (MVDR convention)

    tx = data(:, 2:end).';


    % STEP 5: ARRAY GEOMETRY SETUP %
    % Mic ordering MUST match the OUTPUT script: descending y (most positive
    % first), so that channel index in the CSV maps to the same physical mic
    % position used by the steering vector. Reversed ordering would produce
    % a peak mirror-reflected through the array centre

    n_arr = linspace((Nr-1)/2, -(Nr-1)/2, Nr);
    mic_positions = [x_a * ones(1, Nr); y_a + d_y * n_arr; zeros(1, Nr)];


    % STEP 6: CONVERT TO FREQUENCY DOMAIN %

    fprintf('Converting to frequency domain...\n');

    d_t = time(2) - time(1);
    F_s = 1 / d_t;


    % STEP 7: CHOOSE FFT LENGTH AND ANALYSIS BINS %
    % Pick size_fft so that f_test lands exactly on a bin centre at the
    % target bin spacing — this keeps signal energy in a small number of
    % adjacent bins rather than smearing it across the spectrum, and avoids
    % integrating noise-only bins into the broadband MVDR sum.
    %
    % Bin spacing is F_s / size_fft. Setting size_fft = round(F_s / target)
    % to the nearest integer such that f_test * size_fft / F_s is integer
    % gives an effective bin width close to target_binwidth_hz.

    fprintf('Selecting analysis bins...\n');

    % size_fft chosen so f_test sits on a bin boundary:
    % bin_index_test = round(f_test * size_fft / F_s) is exact integer
    cycles_per_window = round(batch_test_freq / target_binwidth_hz);
    size_fft = round(cycles_per_window * F_s / batch_test_freq);
    actual_binwidth = F_s / size_fft;

    % Keep snapshots reasonably long but well below total recording length
    if size_fft > N / 4
        size_fft = 2^floor(log2(N / 4));
        actual_binwidth = F_s / size_fft;
        fprintf('  Capped size_fft at %d (recording too short for target binwidth)\n', size_fft);
    end

    fft_vec = F_s * (0:(size_fft-1)) / size_fft;

    % Bin indices within ±analysis_halfwidth_hz of f_test
    bin_mask = (fft_vec >= f_band_lo) & (fft_vec <= f_band_hi);
    bin_index = find(bin_mask);
    bin_freqs = fft_vec(bin_index);
    num_bins = length(bin_index);

    if num_bins == 0
        % Fallback — analysis band narrower than one bin width; use single
        % nearest bin to f_test
        [~, nearest] = min(abs(fft_vec(1:floor(size_fft/2)) - batch_test_freq));
        bin_index = nearest;
        bin_freqs = fft_vec(nearest);
        num_bins = 1;
        fprintf('  Analysis band narrower than bin width — using single nearest bin\n');
    end

    fprintf('  size_fft = %d | bin width = %.3f Hz | %d bins (%.1f-%.1f Hz)\n', ...
        size_fft, actual_binwidth, num_bins, min(bin_freqs), max(bin_freqs));


    % STEP 8: CREATE SNAPSHOTS AND CROSS-SPECTRAL MATRIX %

    fprintf('Creating snapshots and cross-spectral matrix...\n');

    window = hanning(size_fft)';
    snapshots = make_snapshots(tx, size_fft, overlap, window);
    maxnum_snap = size(snapshots, 3);

    fprintf('  Created %d snapshots\n', maxnum_snap);

    r = create_csm(snapshots, bin_index, Nr, num_bins, maxnum_snap);
    fprintf('  Created %dx%dx%d cross-spectral matrix\n', size(r));


    % STEP 9: DEFINE GRID SEARCH AREA %

    fprintf('Defining search area...\n');

    system_centre_y = -0.16; %geometric midpoint of full aperture (m)
    [x_scan, y_scan, X_grid, Y_grid, candidate_points, grid_res] = ...
        build_search_grid(batch_test_freq, c_0, source_distance_m, ...
        system_centre_y, grid_pts_per_lambda, margin_fixed, grid_resolution_fixed);

    fprintf('  Search grid: %d x %d points, resolution: %.4f m\n', ...
        length(x_scan), length(y_scan), grid_res);


    % STEP 10: PERFORM MVDR BEAMFORMING %

    fprintf('\n<strong>Running MVDR Beamforming...</strong>\n');
    fprintf('  Regularisation parameter: %.2e\n', loading);

    source_positions = [source_x, source_y, 0];
    source_frequencies = batch_test_freq;

    % 1D scan along y-axis at source x for verification
    fprintf('\n1D Scan (x = %.3f m):\n', source_x);
    plot_1dscan_mvdr(r, mic_positions, source_x, source_y, bin_freqs, c_0, loading);
    exportgraphics(gcf, fullfile(results_folder, ...
        sprintf('%02d_%.0fhz_1dscan.png', batch_test_idx, batch_test_freq)), ...
        'ContentType', 'image');

    % 2D scan of full search area
    fprintf('\n2D Scan:\n');
    [est_x, est_y] = plot_2dscan_mvdr(r, mic_positions, candidate_points, ...
        bin_freqs, c_0, y_scan, x_scan, source_x, source_y, X_grid, Y_grid, loading, total_aperture);
    exportgraphics(gcf, fullfile(results_folder, ...
        sprintf('%02d_%.0fhz_2dscan.png', batch_test_idx, batch_test_freq)), ...
        'ContentType', 'image');
    exportgraphics(gcf, fullfile(results_folder, ...
        sprintf('%02d_%.0fhz_2dscan.pdf', batch_test_idx, batch_test_freq)), ...
        'ContentType', 'vector');


    % STEP 11: PERFORMANCE METRICS %

    fprintf('\n<strong>PERFORMANCE METRICS</strong>\n');

    est_pos = [est_x, est_y];
    true_pos = [source_x, source_y];
    errors = true_pos - est_pos;

    radial_error = norm(errors);
    MSE = sum(errors.^2);
    MSE_percent = 100 * MSE / (grid_res^2);

    true_angle = atan2(source_y - y_a, source_x - x_a);
    est_angle = atan2(est_y - y_a, est_x - x_a);
    angular_error_deg = rad2deg(abs(true_angle - est_angle));
    if angular_error_deg > 180
        angular_error_deg = 360 - angular_error_deg;
    end

    batch_res_radial_error(batch_test_idx) = radial_error;
    batch_res_angular_error(batch_test_idx) = angular_error_deg;

    fprintf('  Estimated position: (%.4f, %.4f) m\n', est_x, est_y);
    fprintf('  Position error: dx = %.4f m, dy = %.4f m\n', errors(1), errors(2));
    fprintf('  Radial error: %.4f m\n', radial_error);
    fprintf('  Mean squared error: %.4f m^2\n', MSE);
    fprintf('  MSE (normalised): %.3f%% of grid resolution\n', MSE_percent);
    fprintf('  Angular error: %.4f deg\n', angular_error_deg);


    % STEP 12: BEAM PATTERN ANALYSIS %

    fprintf('\n<strong>GENERATING BEAM PATTERNS</strong>\n');

    array_centre_2d = [x_a; y_a];
    radius = norm([source_x - x_a, source_y - y_a]);
    fprintf('Using radius = %.3f m (source distance)\n', radius);

    [~, ~, fig_handle, beam_metrics] = beam_pattern_ula(...
        r, mic_positions, array_centre_2d, source_positions(:,1:2), ...
        bin_freqs, c_0, loading, radius, ...
        sprintf('ULA Beam Pattern (f = %.0f Hz)', batch_test_freq));

    exportgraphics(fig_handle, fullfile(results_folder, ...
        sprintf('%02d_%.0fhz_beam.png', batch_test_idx, batch_test_freq)), ...
        'ContentType', 'image');
    close(fig_handle);

    batch_res_beamwidth(batch_test_idx) = beam_metrics.beamwidth_3dB;
    batch_res_sidelobe(batch_test_idx) = beam_metrics.sidelobe_level;
    batch_res_DI(batch_test_idx) = beam_metrics.directivity_index;


    % STEP 13: 1D RAY RESPONSE ANALYSIS %

    % Ray points from array centre through the estimated source position —
    % ray angle is derived from the 2D scan result rather than preset, so
    % the response is sampled along the direction the beamformer actually
    % identifies as the source bearing
    fprintf('\n<strong>1D RAY RESPONSE ANALYSIS</strong>\n');

    ray_angle_deg_est = rad2deg(atan2(est_y - y_a, est_x - x_a));
    fprintf('  Ray angle (through estimate): %.2f deg\n', ray_angle_deg_est);

    [ray_distances, ray_response_db, peak_r] = ray_response_1d(...
        r, mic_positions, array_centre_2d, ray_angle_deg_est, ...
        ray_range_lambdas, ray_num_points, bin_freqs, c_0, loading, lambda);

    batch_res_ray_peak(batch_test_idx) = peak_r;

    fprintf('  True source distance: %.4f m (%.2f lambda)\n', radius, radius / lambda);
    fprintf('  Peak response distance: %.4f m (%.2f lambda)\n', peak_r, peak_r / lambda);

    fig_ray = plot_ray_response(ray_distances, ray_response_db, ...
        radius, lambda, ray_angle_deg_est, batch_test_freq);
    exportgraphics(fig_ray, fullfile(results_folder, ...
        sprintf('%02d_%.0fhz_ray_response.png', batch_test_idx, batch_test_freq)), ...
        'ContentType', 'image');
    exportgraphics(fig_ray, fullfile(results_folder, ...
        sprintf('%02d_%.0fhz_ray_response.pdf', batch_test_idx, batch_test_freq)), ...
        'ContentType', 'vector');
    close(fig_ray);


    fprintf('\n<strong>Test %d Summary (f = %.0f Hz):</strong>\n', ...
        batch_test_idx, batch_test_freq);
    fprintf('  Source: (%.3f, %.3f) m | lambda = %.4f m\n', source_x, source_y, lambda);
    fprintf('  Radial error: %.4f m | Angular error: %.3f deg\n', ...
        radial_error, angular_error_deg);
    fprintf('  Beamwidth: %.1f deg | Sidelobe: %.2f dB | DI: %.2f dB\n', ...
        beam_metrics.beamwidth_3dB, beam_metrics.sidelobe_level, ...
        beam_metrics.directivity_index);
    fprintf('  Ray peak: %.4f m (%.2f lambda)\n', peak_r, peak_r / lambda);

    close all;

end


%% SAVE RESULTS TO CSV %%

fprintf('\n========================================\n');
fprintf('<strong>SAVING RESULTS</strong>\n');
fprintf('========================================\n');

results_table = table(...
    batch_res_frequency_hz, batch_res_lambda_m, ...
    batch_res_source_x, batch_res_source_y, batch_res_grid_resolution, ...
    batch_res_radial_error, batch_res_angular_error, ...
    batch_res_beamwidth, batch_res_sidelobe, batch_res_DI, ...
    batch_res_ray_peak, ...
    'VariableNames', {'Frequency_Hz', 'Lambda_m', ...
    'Source_X_m', 'Source_Y_m', 'Grid_Resolution_m', ...
    'Radial_Error_m', 'Angular_Error_deg', ...
    'Beamwidth_3dB_deg', 'Sidelobe_dB', 'DI_dB', ...
    'Ray_Peak_m'});

fprintf('\nResults table preview:\n');
disp(results_table);

csv_out = fullfile(results_folder, 'freq_sweep_results.csv');
writetable(results_table, csv_out);
fprintf('\nResults saved to: %s\n', csv_out);

mat_out = fullfile(results_folder, 'freq_sweep_results.mat');
save(mat_out, 'results_table', 'batch_test_frequencies', ...
    'source_distance_m', 'source_angle_deg', 'c_0', 'd_y', 'N_mics');
fprintf('Results also saved to: %s\n', mat_out);


%% PUBLICATION FIGURES — FREQUENCY SWEEP RESULTS %%

col_ula = [0.000, 0.447, 0.741]; %blue
lw = 1.8;
ms = 7;
freq_ticks = batch_test_frequencies;


% FIGURE: RADIAL ERROR VS FREQUENCY %

figure('Color', 'w', 'Position', [100 100 520 360]);
semilogy(batch_res_frequency_hz, batch_res_radial_error, '-o', ...
    'Color', col_ula, 'LineWidth', lw, 'MarkerSize', ms, ...
    'MarkerFaceColor', col_ula, 'DisplayName', 'ULA');
xlabel('Frequency (Hz)', 'FontName', 'Times New Roman', 'FontSize', 14);
ylabel('Radial error (m)', 'FontName', 'Times New Roman', 'FontSize', 14);
legend('Location', 'best', 'FontName', 'Times New Roman', 'FontSize', 11, 'Box', 'on');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12, ...
    'XTick', freq_ticks, 'XTickLabelRotation', 0, ...
    'LineWidth', 1.0, 'Box', 'on', 'YGrid', 'on', 'XGrid', 'off');
exportgraphics(gcf, fullfile(results_folder, 'radial_error_vs_frequency.png'), ...
    'ContentType', 'image');


% FIGURE: ANGULAR ERROR VS FREQUENCY %

figure('Color', 'w', 'Position', [100 100 520 360]);
plot(batch_res_frequency_hz, batch_res_angular_error, '-o', ...
    'Color', col_ula, 'LineWidth', lw, 'MarkerSize', ms, ...
    'MarkerFaceColor', col_ula, 'DisplayName', 'ULA');
xlabel('Frequency (Hz)', 'FontName', 'Times New Roman', 'FontSize', 14);
ylabel('Angular error (deg)', 'FontName', 'Times New Roman', 'FontSize', 14);
legend('Location', 'best', 'FontName', 'Times New Roman', 'FontSize', 11, 'Box', 'on');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12, ...
    'XTick', freq_ticks, 'LineWidth', 1.0, 'Box', 'on', ...
    'YGrid', 'on', 'XGrid', 'off');
exportgraphics(gcf, fullfile(results_folder, 'angular_error_vs_frequency.png'), ...
    'ContentType', 'image');


% FIGURE: BEAMWIDTH AND DIRECTIVITY INDEX (2x1 PANEL) %

figure('Color', 'w', 'Position', [100 100 520 600]);

subplot(2,1,1);
plot(batch_res_frequency_hz, batch_res_beamwidth, '-o', ...
    'Color', col_ula, 'LineWidth', lw, 'MarkerSize', ms, ...
    'MarkerFaceColor', col_ula, 'DisplayName', 'ULA');
ylabel('$-3\,\mathrm{dB}$ beamwidth (deg)', 'Interpreter', 'latex', 'FontSize', 13);
legend('Location', 'best', 'FontName', 'Times New Roman', 'FontSize', 11, 'Box', 'on');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12, ...
    'XTick', freq_ticks, 'XTickLabel', {}, ...
    'LineWidth', 1.0, 'Box', 'on', 'YGrid', 'on', 'XGrid', 'off');

subplot(2,1,2);
plot(batch_res_frequency_hz, batch_res_DI, '-o', ...
    'Color', col_ula, 'LineWidth', lw, 'MarkerSize', ms, ...
    'MarkerFaceColor', col_ula, 'DisplayName', 'ULA');
xlabel('Frequency (Hz)', 'FontName', 'Times New Roman', 'FontSize', 14);
ylabel('Directivity index (dB)', 'FontName', 'Times New Roman', 'FontSize', 14);
legend('Location', 'best', 'FontName', 'Times New Roman', 'FontSize', 11, 'Box', 'on');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12, ...
    'XTick', freq_ticks, 'LineWidth', 1.0, 'Box', 'on', ...
    'YGrid', 'on', 'XGrid', 'off');
exportgraphics(gcf, fullfile(results_folder, 'beamwidth_DI_vs_frequency.png'), ...
    'ContentType', 'image');


%% FINAL SUMMARY %%

fprintf('\n========================================\n');
fprintf('<strong>FREQUENCY SWEEP TEST COMPLETE</strong>\n');
fprintf('========================================\n');

fprintf('\n%-12s %-16s %-16s %-16s\n', 'Freq (Hz)', 'Radial (m)', 'Angular (deg)', 'Ray Peak (m)');
fprintf('%s\n', repmat('-', 1, 62));
for i = 1:batch_num_tests
    fprintf('%-12.0f %-16.4f %-16.3f %-16.4f\n', ...
        batch_res_frequency_hz(i), batch_res_radial_error(i), ...
        batch_res_angular_error(i), batch_res_ray_peak(i));
end
fprintf('%s\n', repmat('-', 1, 62));
fprintf('\nAll results saved to: %s\n\n', results_folder);


%% FUNCTION DEFINITIONS %%

% FUNCTION: MAKE_SNAPSHOTS
% Make time-domain snapshots of input signal
function snapshots = make_snapshots(tx, size_fft, overlap, window)

    step = round(overlap * size_fft);
    num_snap = floor((size(tx, 2) - size_fft) / step) + 1;

    start_idx = 1 + (0:(num_snap-1)) * step;
    idx_matrix = start_idx + (0:size_fft-1)';

    snapshots = tx(:, idx_matrix);
    snapshots = reshape(snapshots, size(tx, 1), size_fft, num_snap);
    snapshots = snapshots .* reshape(window, [1, size_fft, 1]);
end


% FUNCTION: CREATE_CSM
% Convert each snapshot block to frequency domain and compute cross-spectral matrix
function r = create_csm(snapshots, bin_index, Nr, num_bins, maxnum_snap)

    fx = fft(snapshots, [], 2); %(Nr x size_fft x num_snap)
    fx = fx(:, bin_index, :); %(Nr x num_bins x num_snap)

    r = zeros(Nr, Nr, num_bins);

    for jf = 1:num_bins
        fx1 = squeeze(fx(:, jf, :)); %(Nr x num_snap)
        r(:, :, jf) = fx1 * fx1';
    end

    r = r / maxnum_snap;
end


% FUNCTION: MVDR_BEAMFORMING
% Eigenvalue decomposition with adaptive CBF-based diagonal loading
function response_db = mvdr_beamforming(r, mic_positions, candidate_points, ...
    bin_freqs, c_0, loading)

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

            % Adaptive regularisation based on CBF output
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
    bin_freqs, c_0, y_scan, x_scan, source_x, source_y, X_grid, Y_grid, loading, total_aperture)

    response_db = mvdr_beamforming(r, mic_positions, candidate_points, ...
        bin_freqs, c_0, loading);
    grid_response = reshape(response_db, length(y_scan), length(x_scan));

    figure('Color', 'w', 'Position', [100 100 520 440]);
    imagesc(x_scan, y_scan, grid_response);
    axis xy;
    axis tight;
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
    grid off;
    hold on;

    [est_x, est_y, max_val] = refine_peak_2d(grid_response, x_scan, y_scan);

    fprintf('  Max response: %.2f dB at (%.3f, %.3f) m\n', max_val, est_x, est_y);

    % True source — black cross
    plot(source_x, source_y, 'kx', ...
        'MarkerSize', 5, 'LineWidth', 1.5, 'DisplayName', 'True source');

    % Estimated source — black hollow circle
    plot(est_x, est_y, 'ko', ...
        'MarkerSize', 10, 'LineWidth', 1.5, ...
        'MarkerFaceColor', 'none', 'DisplayName', 'Estimated source');

    % Array aperture — vertical line spanning full array extent with endpoint ticks
    % represents the physical ULA along the y-axis at x = x_a
    y_aperture_top = max(mic_positions(2,:));
    y_aperture_bot = max(mic_positions(2,:)) - total_aperture;
    x_array = mean(mic_positions(1,:));
    tick_half = 0.03; %half-width of endpoint ticks (m)
    col_array = [0.7 0.7 0.7]; %medium grey

    plot([x_array, x_array], [y_aperture_bot, y_aperture_top], '-', ...
        'LineWidth', 1.5, 'Color', col_array, 'DisplayName', 'Array');
    plot([x_array - tick_half, x_array + tick_half], ...
        [y_aperture_top, y_aperture_top], '-', ...
        'LineWidth', 1.5, 'Color', col_array, 'HandleVisibility', 'off');
    plot([x_array - tick_half, x_array + tick_half], ...
        [y_aperture_bot, y_aperture_bot], '-', ...
        'LineWidth', 1.5, 'Color', col_array, 'HandleVisibility', 'off');

    legend('Location', 'northeast', 'FontName', 'Times New Roman', ...
        'FontSize', 10, 'Color', 'w', 'TextColor', 'k', 'Box', 'on');
end


% FUNCTION: PLOT_1DSCAN_MVDR
% 1D MVDR scan along y-axis at fixed x equal to source x coordinate
function plot_1dscan_mvdr(r, mic_positions, source_x, source_y, bin_freqs, c_0, loading)

    y_scan_line = linspace(source_y - 0.5, source_y + 0.5, 200);
    z_fixed = 0;
    candidate_points = [source_x * ones(length(y_scan_line), 1), ...
                        y_scan_line(:), ...
                        z_fixed * ones(length(y_scan_line), 1)];

    response_line = mvdr_beamforming(r, mic_positions, candidate_points, ...
        bin_freqs, c_0, loading);

    [max_response, max_idx] = max(response_line);
    y_est = y_scan_line(max_idx);

    fprintf('  Max response: %.2f dB at y = %.3f m\n', max_response, y_est);

    figure('Name', '1D MVDR Scan');
    plot(y_scan_line, response_line, 'LineWidth', 2);
    xlabel('$y$ position (m)', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('Beamformer response (dB)', 'FontName', 'Times New Roman', 'FontSize', 14);
    grid on;
    hold on;
    xline(source_y, 'r--', 'LineWidth', 1.5, 'DisplayName', 'True source');
    legend('Beamformer response', 'True source', 'FontName', 'Times New Roman');
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
end


% FUNCTION: RAY_RESPONSE_1D
% Compute MVDR response along a 1D radial ray from the array centre at a
% fixed bearing. Returns the normalised response in dB and the sub-sample
% peak distance via parabolic interpolation
function [ray_distances, ray_response_db, peak_r] = ray_response_1d(...
    r, mic_positions, array_centre_2d, ray_angle_deg, ray_range_lambdas, ...
    ray_num_points, bin_freqs, c_0, loading, lambda)

    r_min = ray_range_lambdas(1) * lambda;
    r_max = ray_range_lambdas(2) * lambda;
    ray_distances = linspace(r_min, r_max, ray_num_points);

    % Build candidate points along the ray — originating at the array centre
    % and sampled radially outward at the specified bearing
    ray_x = array_centre_2d(1) + ray_distances * cosd(ray_angle_deg);
    ray_y = array_centre_2d(2) + ray_distances * sind(ray_angle_deg);
    ray_points = [ray_x(:), ray_y(:), zeros(ray_num_points, 1)];

    ray_response_db = mvdr_beamforming(r, mic_positions, ray_points, ...
        bin_freqs, c_0, loading);

    peak_r = refine_peak_1d(ray_distances, ray_response_db);
end


% FUNCTION: REFINE_PEAK_1D
% Sub-sample peak localisation on a 1D response via parabolic interpolation
function peak_r = refine_peak_1d(distances, response_db)

    [~, pk_idx] = max(response_db);

    if pk_idx < 2 || pk_idx > length(distances) - 1
        peak_r = distances(pk_idx);
        return;
    end

    dr = distances(2) - distances(1);
    f_m = response_db(pk_idx - 1);
    f_0 = response_db(pk_idx);
    f_p = response_db(pk_idx + 1);

    denom = f_m - 2 * f_0 + f_p;
    if abs(denom) > eps
        delta_i = (f_m - f_p) / (2 * denom);
    else
        delta_i = 0;
    end

    delta_i = max(-0.5, min(0.5, delta_i));
    peak_r = distances(pk_idx) + delta_i * dr;
end


% FUNCTION: PLOT_RAY_RESPONSE
% Plot normalised MVDR response along 1D ray with true distance marker
function fig_handle = plot_ray_response(ray_distances, ray_response_db, ...
    source_distance_m, lambda, ray_angle_deg, test_freq)

    col_ula = [0.000, 0.447, 0.741]; %blue — matches ULA plotting style
    fn = 'Times New Roman';

    ray_distances_lambda = ray_distances / lambda;

    fig_handle = figure('Color', 'w', 'Position', [100 100 520 360]);
    plot(ray_distances_lambda, ray_response_db, '-', ...
        'Color', col_ula, 'LineWidth', 1.8, ...
        'DisplayName', 'ULA');
    hold on;
    xline(source_distance_m / lambda, 'k--', 'LineWidth', 1.2, ...
        'DisplayName', 'True distance');
    hold off;

    xlabel('Distance along ray ($r / \lambda$)', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('Normalised response (dB)', 'FontName', fn, 'FontSize', 14);
    title(sprintf('Ray response at %.1f deg (f = %.0f Hz)', ray_angle_deg, test_freq), ...
        'FontName', fn, 'FontSize', 13);
    legend('Location', 'southeast', 'FontName', fn, 'FontSize', 11, 'Box', 'on');
    set(gca, 'FontName', fn, 'FontSize', 12, ...
        'LineWidth', 1.0, 'Box', 'on', 'YGrid', 'on', 'XGrid', 'off');
end


% FUNCTION: BEAM_PATTERN_ULA
% Compute near-field beam pattern for ULA MVDR
% Evaluates MVDR response at each angular direction at fixed range
% Returns beam_metrics struct for batch results storage
function [theta_deg, beam_pattern, fig_handle, beam_metrics] = beam_pattern_ula(...
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

    % Polar plot
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
            source_vec = source_positions(src, :)' - array_centre;
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

    % Beam pattern statistics
    fprintf('\n<strong>Beam Pattern Statistics:</strong>\n');
    [~, peak_idx] = max(beam_pattern);
    fprintf('  Peak response: 0.00 dB at %.1f degrees\n', theta_deg(peak_idx));

    beam_metrics = struct();

    % -3 dB beamwidth — identify the contiguous region containing the peak
    threshold_3db = max(beam_pattern) - 3;
    above = beam_pattern(:) >= threshold_3db;
    d_above = diff([0; above; 0]);
    starts = find(d_above == 1);
    ends_arr = find(d_above == -1) - 1;

    beam_metrics.beamwidth_3dB = NaN;
    main_lobe_region = find(starts <= peak_idx & ends_arr >= peak_idx, 1);
    if ~isempty(main_lobe_region)
        bw = theta_deg(ends_arr(main_lobe_region)) - theta_deg(starts(main_lobe_region));
        if bw < 0
            bw = bw + 360;
        end
        beam_metrics.beamwidth_3dB = bw;
        fprintf('  -3 dB beamwidth: %.1f degrees (%.1f deg to %.1f deg)\n', bw, ...
            theta_deg(starts(main_lobe_region)), theta_deg(ends_arr(main_lobe_region)));
    else
        fprintf('  WARNING: could not identify main lobe region\n');
    end

    % Maximum sidelobe level
    beam_metrics.sidelobe_level = NaN;
    bp_col = beam_pattern(:);
    if any(~above)
        beam_metrics.sidelobe_level = max(bp_col(~above));
        fprintf('  Maximum sidelobe level: %.2f dB\n', beam_metrics.sidelobe_level);
    end

    % Directivity index
    beam_linear = 10.^(beam_pattern(:)/10);
    integral_val = trapz(theta_rad(:), beam_linear);
    DI = 10 * log10(2 * pi / integral_val);
    beam_metrics.directivity_index = DI;
    fprintf('  Directivity index: %.2f dB\n', DI);

    fprintf('\n');
end


% FUNCTION: BUILD_SEARCH_GRID
% Construct a consistent 2D search grid for MVDR beamforming
% Written by L Marshall — grid consistency fix for InterNoise 2026 scripts
function [x_scan, y_scan, X_grid, Y_grid, candidate_points, grid_res] = ...
    build_search_grid(test_freq, c_0, source_distance, array_centre_y, ...
    grid_pts_per_lambda, margin_fixed, grid_resolution_fixed)

    % Grid step — fixed override takes priority over wavelength-scaled value
    if isempty(grid_resolution_fixed)
        lambda = c_0 / test_freq;
        grid_resolution = lambda / grid_pts_per_lambda;
    else
        grid_resolution = grid_resolution_fixed;
    end

    max_extent = source_distance + margin_fixed;
    x_half_extent = max_extent;
    y_half_extent = max_extent;

    x_scan_points = round(2 * x_half_extent / grid_resolution) + 1;
    y_scan_points = round(2 * y_half_extent / grid_resolution) + 1;

    % Centre on system geometric midpoint
    x_centre = 0;
    y_centre = array_centre_y;

    x_scan = linspace(x_centre - x_half_extent, x_centre + x_half_extent, x_scan_points);
    y_scan = linspace(y_centre - y_half_extent, y_centre + y_half_extent, y_scan_points);

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

    denom_x = fx_m - 2*fx_0 + fx_p;
    if abs(denom_x) > eps
        delta_ix = (fx_m - fx_p) / (2 * denom_x);
    else
        delta_ix = 0; %flat region — no refinement possible
    end

    % Step 4: parabolic interpolation along y at x = ix_pk
    fy_m = grid_response(iy_pk - 1, ix_pk);
    fy_0 = grid_response(iy_pk, ix_pk);
    fy_p = grid_response(iy_pk + 1, ix_pk);

    denom_y = fy_m - 2*fy_0 + fy_p;
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