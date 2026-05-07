clc; clear all; close all;

% Batch frequency sweep testing — Acoustic Vector Array MVDR
% Runs signal generation and AVA MVDR processing at a fixed source distance
% across a range of test frequencies for localisation performance comparison
%
% Sweep frequencies: 630, 800, 1000, 1250, 1600, 2000, 2500 Hz
% Fixed source distance: 0.4 m at 180 deg from +X axis
%
% For each frequency this script:
%   (1) Generates a synthetic signal via SUPER_VA_OUTPUT_BATCH.m
%   (2) Runs MVDR beamforming and superimposed localisation
%   (3) Extracts radial error, angular error, and beam pattern metrics
%   (4) Saves all results to CSV and .mat for post-processing
%
% NOTE: The superimposed beam pattern is only valid in simulation where both
% arrays observe the same stationary source field simultaneously. Per-angle
% linear power responses from each array are multiplied before normalisation,
% yielding a narrower main lobe than either sub-array alone. Not valid for
% sequentially recorded experimental data.
%
% Written by L Marshall 29/03/2026
% Corrected FFT bin selection by L Marshall 29/04/2026
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
%   - grid_resolution_fixed parameter added — set to a positive value to
%     override the frequency-dependent λ/N grid with a constant step,
%     useful for isolating grid quantisation effects from genuine
%     frequency-dependent performance variation


%% BATCH TEST PARAMETERS %%

% Test frequencies (Hz)
batch_test_frequencies = [630 800 1000 1250 1600 2000 2500 3150 4000];
batch_num_tests = length(batch_test_frequencies);

% Fixed source geometry
source_distance_m = 0.4; %fixed source distance (m)
source_angle_deg = 180; %source angle from +X axis (deg)
c_0 = 340; %speed of sound (m/s)
rho_0 = 1.02; %air density at STP

% Array parameters
delta_fixed = 0.04; %MEMS colocation spacing (m)
N_a = 1; %number of independent arrays
N_v = 1; %vector sensors per array
d_y = 0.1; %vector sensor y-axis spacing (m)

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
grid_resolution_fixed = [0.01]; %override grid step (m); [] = use λ/grid_pts_per_lambda
margin_fixed = 0.5; %fixed physical search margin (m)

% Create results folder
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
results_folder = sprintf('freq_sweep_AVA_%s', timestamp);
mkdir(results_folder);

fprintf('\n<strong>BATCH FREQUENCY SWEEP TEST - AVA MVDR</strong>\n');
fprintf('Testing %d frequencies: ', batch_num_tests);
fprintf('%.0f  ', batch_test_frequencies);
fprintf('Hz\n');
fprintf('Fixed source distance: %.3f m at %.0f deg\n', source_distance_m, source_angle_deg);
fprintf('Array 2 separation: fixed at (0, -0.32) m\n');
fprintf('Fixed delta: %.4f m\n', delta_fixed);
if isempty(grid_resolution_fixed)
    fprintf('Grid: %d pts/lambda (frequency-dependent) | Margin: %.2f m\n', ...
        grid_pts_per_lambda, margin_fixed);
else
    fprintf('Grid: %.4f m fixed | Margin: %.2f m\n', ...
        grid_resolution_fixed, margin_fixed);
end
fprintf('Results folder: %s\n\n', results_folder);


%% INITIALISE RESULTS STORAGE %%

batch_res_frequency_hz = zeros(batch_num_tests, 1);
batch_res_lambda_m = zeros(batch_num_tests, 1);
batch_res_source_x = zeros(batch_num_tests, 1);
batch_res_source_y = zeros(batch_num_tests, 1);
batch_res_grid_resolution = zeros(batch_num_tests, 1);

% Individual array localisation metrics
batch_res_radial_error_a1 = zeros(batch_num_tests, 1);
batch_res_angular_error_a1 = zeros(batch_num_tests, 1);
batch_res_radial_error_a2 = zeros(batch_num_tests, 1);
batch_res_angular_error_a2 = zeros(batch_num_tests, 1);

% Superimposed localisation metrics
batch_res_radial_error_super = zeros(batch_num_tests, 1);
batch_res_angular_error_super = zeros(batch_num_tests, 1);

% Individual array beam pattern metrics
batch_res_beamwidth_a1 = zeros(batch_num_tests, 1);
batch_res_beamwidth_a2 = zeros(batch_num_tests, 1);
batch_res_sidelobe_a1 = zeros(batch_num_tests, 1);
batch_res_sidelobe_a2 = zeros(batch_num_tests, 1);
batch_res_DI_a1 = zeros(batch_num_tests, 1);
batch_res_DI_a2 = zeros(batch_num_tests, 1);

% Superimposed beam pattern metrics — simulation only
batch_res_beamwidth_super = zeros(batch_num_tests, 1);
batch_res_sidelobe_super = zeros(batch_num_tests, 1);
batch_res_DI_super = zeros(batch_num_tests, 1);


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

    % Array 2 position — fixed physical separation
    x_a2 = 0;
    y_a2 = -0.32;

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
    fprintf('Array 2 centre: (%.3f, %.3f) m | Separation: %.3f m (%.2f lambda)\n', ...
        x_a2, y_a2, norm([x_a2, y_a2]), norm([x_a2, y_a2]) / lambda);
    fprintf('Analysis band: %.0f-%.0f Hz\n', f_band_lo, f_band_hi);


    % STEP 1: CLEAR PREVIOUS RUN VARIABLES %

    vars_to_clear = {'source_positions', 'source_frequencies', 'source_amplitudes', ...
        'vs_centres_arrays', 'r_vs_arrays', 'response_db_arrays', ...
        'response_db_combined', 'response_db_norm_combined', ...
        'all_cases_results', 'all_results'};
    for v = 1:length(vars_to_clear)
        if exist(vars_to_clear{v}, 'var')
            clear(vars_to_clear{v});
        end
    end


    % STEP 2: RUN SIGNAL GENERATION SCRIPT %

    fprintf('\nGenerating signal at r=%.3f m, f=%.0f Hz...\n', ...
        source_distance_m, batch_test_freq);

    batch_test_delta = delta_fixed;
    batch_test_source_x = source_x;
    batch_test_source_y = source_y;
    batch_csv_name = fullfile(results_folder, ...
        sprintf('signal_%02d_%.0fhz.csv', batch_test_idx, batch_test_freq));

    fprintf('  batch_test_freq = %.0f Hz\n', batch_test_freq);
    fprintf('  batch_test_delta = %.4f m\n', batch_test_delta);
    fprintf('  batch_test_source_x = %.3f m\n', batch_test_source_x);
    fprintf('  batch_test_source_y = %.3f m\n', batch_test_source_y);
    fprintf('  batch_csv_name = %s\n', batch_csv_name);

    run('SUPER_VA_OUTPUT_BATCH.m');
    fprintf('  Signal generation complete\n');


    % STEP 3: LOAD SIGNAL FROM CSV %
    % CSV format: [time, mic_1, mic_2, ..., mic_N_m] — single signal per channel
    % Note: SUPER_VA_OUTPUT generates with N_a = 1 internally, so the CSV
    % from the OUTPUT script contains only 4 mics (one AVS). For the
    % two-array configuration here we re-generate the second array's
    % signal by reading per-array channels from the CSV; the OUTPUT script
    % needs N_a = 2 set internally for full dual-array data.

    fprintf('\nLoading signal from: %s\n', csv_filename);

    data_table = readtable(csv_filename);
    data = table2array(data_table);
    time = data(:, 1);
    N = size(data, 1);
    N_ma = N_v * 4; %microphones per array
    N_m_total = size(data, 2) - 1;

    if N_m_total ~= N_a * N_ma
        warning('CSV holds %d mics; expected %d for N_a=%d, N_v=%d. Check OUTPUT script N_a setting.', ...
            N_m_total, N_a * N_ma, N_a, N_v);
    end

    fprintf('  Loaded %d samples | %d arrays | %d mics each\n', N, N_a, N_ma);


    % STEP 4: ASSEMBLE PER-ARRAY SIGNAL MATRICES %
    % Rows = sensors, columns = time samples (MVDR convention)

    tx_vs_arrays = cell(N_a, 1);
    for array = 1:N_a
        col_start = 1 + (array - 1) * N_ma + 1;
        col_end = col_start + N_ma - 1;
        tx_vs_arrays{array} = data(:, col_start:col_end).';
    end


    % STEP 5: ARRAY GEOMETRY SETUP %

    A1_coord = [0; 0; 0];
    A2_coord = [x_a2; y_a2; 0];
    array_centres = [A1_coord, A2_coord];

    mic_offsets = delta_fixed / 2 * [-1, 1, 1, -1; -1, -1, 1, 1; 0, 0, 0, 0];

    vs_centres_arrays = cell(N_a, 1);
    mic_positions = zeros(3, N_a * N_ma);

    for array = 1:N_a
        vs_centres = zeros(3, N_v);
        for vs = 1:N_v
            vs_centres(:, vs) = array_centres(:, array) + [0; (vs - 1) * d_y; 0];
            idx = (array - 1) * N_ma + (vs - 1) * 4 + (1:4);
            mic_positions(:, idx) = vs_centres(:, vs) + mic_offsets;
        end
        vs_centres_arrays{array} = vs_centres;
    end

    fprintf('  Array separation: %.3f m\n', norm(diff(array_centres, 1, 2)));


    % STEP 6: CONVERT TO FREQUENCY DOMAIN %

    fprintf('Converting to frequency domain...\n');

    d_t = time(2) - time(1);
    F_s = 1 / d_t;


    % STEP 7: CHOOSE FFT LENGTH AND ANALYSIS BINS %
    % Pick size_fft so that f_test lands exactly on a bin centre at the
    % target bin spacing — keeps signal energy in a small number of
    % adjacent bins rather than smearing it across the spectrum, and
    % avoids integrating noise-only bins into the broadband MVDR sum.

    fprintf('Selecting analysis bins...\n');

    cycles_per_window = round(batch_test_freq / target_binwidth_hz);
    size_fft = round(cycles_per_window * F_s / batch_test_freq);
    actual_binwidth = F_s / size_fft;

    % Cap at recording length for sufficient snapshot count
    if size_fft > N / 4
        size_fft = 2^floor(log2(N / 4));
        actual_binwidth = F_s / size_fft;
        fprintf('  Capped size_fft at %d (recording too short for target binwidth)\n', size_fft);
    end

    fft_vec = F_s * (0:(size_fft - 1)) / size_fft;

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


    % STEP 8: CREATE AVS SNAPSHOTS AND CROSS-SPECTRAL MATRIX %

    fprintf('Creating snapshots and cross-spectral matrix...\n');

    window = hanning(size_fft)';
    snapshots_vs_arrays = cell(N_a, 1);
    r_vs_arrays = cell(N_a, 1);

    for array = 1:N_a
        snapshots_vs_arrays{array} = make_snapshots(tx_vs_arrays{array}, ...
            size_fft, overlap, window);
        r_vs_arrays{array} = create_vs_csm(snapshots_vs_arrays{array}, bin_index, ...
            delta_fixed, bin_freqs, rho_0, N_v, num_bins, c_0);
        fprintf('  Array %d: created %d snapshots | CSM %dx%dx%d\n', array, ...
            size(snapshots_vs_arrays{array}, 3), size(r_vs_arrays{array}));
    end


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
    num_sources = 1;

    response_db_arrays = cell(N_a, 1);
    for array = 1:N_a
        fprintf('  Processing Array %d...\n', array);
        response_db_arrays{array} = mvdr_vs_beamforming(r_vs_arrays{array}, ...
            vs_centres_arrays{array}, candidate_points, bin_freqs, c_0, rho_0, ...
            loading, num_bins, N_v);
        fprintf('    Array %d response range: [%.2f, %.2f] dB\n', array, ...
            min(response_db_arrays{array}), max(response_db_arrays{array}));
    end

    % 2D scans for individual arrays
    fprintf('\n2D Scans - Individual Arrays:\n');
    for array = 1:N_a
        fprintf('  Plotting Array %d response...\n', array);
        plot_2dscan_mvdr(r_vs_arrays{array}, vs_centres_arrays{array}, ...
            candidate_points, bin_freqs, c_0, rho_0, y_scan, x_scan, ...
            X_grid, Y_grid, loading, num_bins, N_v, num_sources, source_positions, ...
            response_db_arrays{array}, array, grid_res, vs_centres_arrays, N_a);
        exportgraphics(gcf, fullfile(results_folder, ...
            sprintf('%02d_%.0fhz_array%d_2dscan.png', batch_test_idx, batch_test_freq, array)), ...
            'ContentType', 'image');
        exportgraphics(gcf, fullfile(results_folder, ...
            sprintf('%02d_%.0fhz_array%d_2dscan.pdf', batch_test_idx, batch_test_freq, array)), ...
            'ContentType', 'vector');
    end


    % STEP 11: COMBINE AND SUPERIMPOSE ARRAY RESPONSES %

    fprintf('\nCombining array responses...\n');

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

    % Mean-normalised combination with ray intersection as source estimate
    [est_positions_combined, ~, response_db_norm_combined] = ...
        superimposed_localisation(response_db_combined, X_grid, Y_grid, y_scan, x_scan, ...
        num_sources, source_positions, N_a, grid_res, response_db_arrays, vs_centres_arrays);

    exportgraphics(gcf, fullfile(results_folder, ...
        sprintf('%02d_%.0fhz_superimposed.png', batch_test_idx, batch_test_freq)), ...
        'ContentType', 'image');
    exportgraphics(gcf, fullfile(results_folder, ...
        sprintf('%02d_%.0fhz_superimposed.pdf', batch_test_idx, batch_test_freq)), ...
        'ContentType', 'vector');


    % STEP 12: PERFORMANCE METRICS %

    fprintf('\n<strong>PERFORMANCE METRICS</strong>\n');

    response_names = [arrayfun(@(i) sprintf('Array %d', i), 1:N_a, 'UniformOutput', false), ...
        'Superimposed'];
    response_db_arrays_norm = cellfun(@(r) r - max(r), response_db_arrays, 'UniformOutput', false);
    responses_to_analyse = [response_db_arrays_norm', {response_db_norm_combined}];

    % Angular reference points: individual array centres + system midpoint
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

    % Store localisation results for this frequency step
    batch_res_radial_error_a1(batch_test_idx) = all_cases_results{1}.radial_error(1);
    batch_res_angular_error_a1(batch_test_idx) = all_cases_results{1}.angular_error(1);
    batch_res_radial_error_a2(batch_test_idx) = all_cases_results{2}.radial_error(1);
    batch_res_angular_error_a2(batch_test_idx) = all_cases_results{2}.angular_error(1);
    batch_res_radial_error_super(batch_test_idx) = all_cases_results{3}.radial_error(1);
    batch_res_angular_error_super(batch_test_idx) = all_cases_results{3}.angular_error(1);

    all_results = all_cases_results{end};


    % STEP 12b: SIMULTANEOUS CSM COMPARISON (SIMULATION ONLY) %
    % Constructs the full 6x6 CSM including cross-array coherence terms.
    % Valid in simulation only — both arrays observe the same field
    % simultaneously. Presented as a performance ceiling comparison.

    % fprintf('\n<strong>Simultaneous CSM comparison (simulation only)...</strong>\n');
    % 
    % r_vs_simultaneous = create_vs_csm_simultaneous(...
    %     snapshots_vs_arrays{1}, snapshots_vs_arrays{2}, ...
    %     bin_index, delta_fixed, bin_freqs, rho_0, N_v, num_bins, c_0);
    % 
    % fprintf('  Full CSM: %dx%dx%d\n', size(r_vs_simultaneous));
    % 
    % response_db_simultaneous = mvdr_vs_beamforming_simultaneous(...
    %     r_vs_simultaneous, vs_centres_arrays{1}, vs_centres_arrays{2}, ...
    %     candidate_points, bin_freqs, c_0, rho_0, loading, num_bins, N_v);
    % 
    % grid_response_sim = reshape(response_db_simultaneous, length(y_scan), length(x_scan));
    % [est_x_sim, est_y_sim, ~] = refine_peak_2d(grid_response_sim, x_scan, y_scan);
    % 
    % figure('Color', 'w', 'Position', [100 100 520 440]);
    % imagesc(x_scan, y_scan, grid_response_sim);
    % axis xy; axis tight;
    % colormap(gca, 'jet');
    % cb = colorbar;
    % cb.Label.String = 'Normalised Response (dB)';
    % cb.Label.Interpreter = 'latex';
    % xlabel('$x$ (m)', 'Interpreter', 'latex', 'FontSize', 14);
    % ylabel('$y$ (m)', 'Interpreter', 'latex', 'FontSize', 14);
    % set(gca, 'FontName', 'Times New Roman', 'FontSize', 12, ...
    %     'LineWidth', 1.0, 'Box', 'on', 'Layer', 'top');
    % grid off; hold on;
    % plot(source_x, source_y, 'kx', 'MarkerSize', 5, 'LineWidth', 1.5, ...
    %     'DisplayName', 'True source');
    % plot(est_x_sim, est_y_sim, 'ko', 'MarkerSize', 10, 'LineWidth', 1.5, ...
    %     'MarkerFaceColor', 'none', 'DisplayName', 'Estimated source');
    % for array = 1:N_a
    %     ac = mean(vs_centres_arrays{array}(1:2,:), 2);
    %     if array == 1
    %         plot(ac(1), ac(2), 'ks', 'MarkerSize', 9, 'LineWidth', 1, ...
    %             'MarkerFaceColor', 'w', 'DisplayName', 'Array');
    %     else
    %         plot(ac(1), ac(2), 'ks', 'MarkerSize', 9, 'LineWidth', 1, ...
    %             'MarkerFaceColor', 'w', 'HandleVisibility', 'off');
    %     end
    % end
    % legend('Location', 'northeast', 'FontName', 'Times New Roman', ...
    %     'FontSize', 10, 'Color', 'w', 'TextColor', 'k', 'Box', 'on');
    % exportgraphics(gcf, fullfile(results_folder, ...
    %     sprintf('%02d_%.0fhz_simultaneous_csm.png', batch_test_idx, batch_test_freq)), ...
    %     'ContentType', 'image');
    % 
    % err_sim = norm([est_x_sim - source_x, est_y_sim - source_y]);
    % fprintf('  Simultaneous CSM: estimated (%.4f, %.4f) m | radial error = %.4f m\n', ...
    %     est_x_sim, est_y_sim, err_sim);
    % fprintf('  Sequential superimposed: radial error = %.4f m\n', ...
    %     batch_res_radial_error_super(batch_test_idx));


    % STEP 13: BEAM PATTERN ANALYSIS %

    fprintf('\n<strong>GENERATING BEAM PATTERNS</strong>\n');

    % Radius = median source distance across all arrays
    distances_all = [];
    for array = 1:N_a
        array_centre_2d = mean(vs_centres_arrays{array}(1:2,:), 2);
        distances = sqrt(sum((source_positions(:,1:2) - array_centre_2d').^2, 2));
        distances_all = [distances_all; distances];
    end
    radius = median(distances_all);
    fprintf('Using radius = %.3f m (median source distance across all arrays)\n', radius);

    % Per-angle linear responses are stored for the superimposed beam pattern
    mvdr_responses_per_array = cell(N_a, 1);

    for array = 1:N_a
        array_centre_2d = mean(vs_centres_arrays{array}(1:2,:), 2);

        [~, ~, fig_handle, beam_metrics, mvdr_resp_linear] = ...
            compute_nearfield_beam_pattern(...
            r_vs_arrays{array}, vs_centres_arrays{array}, array_centre_2d, ...
            source_positions(:,1:2), bin_freqs, c_0, rho_0, loading, num_bins, N_v, ...
            radius, sprintf('Array %d Near-Field Beam Pattern (f = %.0f Hz)', ...
            array, batch_test_freq));

        exportgraphics(fig_handle, fullfile(results_folder, ...
            sprintf('%02d_%.0fhz_array%d_beam.png', batch_test_idx, batch_test_freq, array)), ...
            'ContentType', 'image');
        exportgraphics(fig_handle, fullfile(results_folder, ...
            sprintf('%02d_%.0fhz_array%d_beam.pdf', batch_test_idx, batch_test_freq, array)), ...
            'ContentType', 'vector');
        close(fig_handle);

        mvdr_responses_per_array{array} = mvdr_resp_linear;

        if array == 1
            batch_res_beamwidth_a1(batch_test_idx) = beam_metrics.beamwidth_3dB;
            batch_res_sidelobe_a1(batch_test_idx) = beam_metrics.sidelobe_level;
            batch_res_DI_a1(batch_test_idx) = beam_metrics.directivity_index;
        else
            batch_res_beamwidth_a2(batch_test_idx) = beam_metrics.beamwidth_3dB;
            batch_res_sidelobe_a2(batch_test_idx) = beam_metrics.sidelobe_level;
            batch_res_DI_a2(batch_test_idx) = beam_metrics.directivity_index;
        end
    end


    % STEP 14: SUPERIMPOSED BEAM PATTERN (SIMULATION ONLY) %
    % Per-angle linear power responses from each array are multiplied before
    % normalisation — consistent with the mean-normalised spatial map combination.
    % High response persists only where both arrays agree, yielding a narrower
    % effective main lobe than either sub-array alone.
    % 
    % fprintf('\n<strong>Computing Superimposed Beam Pattern (simulation only)...</strong>\n');
    % 
    % system_centre_2d = mean(array_geometric_centres, 2);
    % 
    % [~, ~, fig_handle_super, beam_metrics_super] = ...
    %     compute_superimposed_beam_pattern(...
    %     mvdr_responses_per_array, system_centre_2d, source_positions(:,1:2), ...
    %     radius, sprintf('Superimposed AVS Beam Pattern (f = %.0f Hz)', batch_test_freq));
    % 
    % exportgraphics(fig_handle_super, fullfile(results_folder, ...
    %     sprintf('%02d_%.0fhz_superimposed_beam.png', batch_test_idx, batch_test_freq)), ...
    %     'ContentType', 'image');
    % close(fig_handle_super);
    % 
    % batch_res_beamwidth_super(batch_test_idx) = beam_metrics_super.beamwidth_3dB;
    % batch_res_sidelobe_super(batch_test_idx) = beam_metrics_super.sidelobe_level;
    % batch_res_DI_super(batch_test_idx) = beam_metrics_super.directivity_index;
    % 
    % fprintf('\n<strong>Test %d Summary (f = %.0f Hz):</strong>\n', ...
    %     batch_test_idx, batch_test_freq);
    % fprintf('  Source: (%.3f, %.3f) m | lambda = %.4f m\n', source_x, source_y, lambda);
    % fprintf('  Radial error: A1 = %.4f m | A2 = %.4f m | Super = %.4f m\n', ...
    %     batch_res_radial_error_a1(batch_test_idx), ...
    %     batch_res_radial_error_a2(batch_test_idx), ...
    %     batch_res_radial_error_super(batch_test_idx));
    % fprintf('  Angular error: A1 = %.3f deg | A2 = %.3f deg | Super = %.3f deg\n', ...
    %     batch_res_angular_error_a1(batch_test_idx), ...
    %     batch_res_angular_error_a2(batch_test_idx), ...
    %     batch_res_angular_error_super(batch_test_idx));
    % fprintf('  Beamwidth: A1 = %.1f deg | A2 = %.1f deg | Super = %.1f deg\n', ...
    %     batch_res_beamwidth_a1(batch_test_idx), ...
    %     batch_res_beamwidth_a2(batch_test_idx), ...
    %     batch_res_beamwidth_super(batch_test_idx));

    close all;

end


%% SAVE RESULTS TO CSV %%

fprintf('\n========================================\n');
fprintf('<strong>SAVING RESULTS</strong>\n');
fprintf('========================================\n');

results_table = table(...
    batch_res_frequency_hz, batch_res_lambda_m, ...
    batch_res_source_x, batch_res_source_y, batch_res_grid_resolution, ...
    batch_res_radial_error_a1, batch_res_angular_error_a1, ...
    batch_res_radial_error_a2, batch_res_angular_error_a2, ...
    batch_res_radial_error_super, batch_res_angular_error_super, ...
    batch_res_beamwidth_a1, batch_res_beamwidth_a2, batch_res_beamwidth_super, ...
    batch_res_sidelobe_a1, batch_res_sidelobe_a2, batch_res_sidelobe_super, ...
    batch_res_DI_a1, batch_res_DI_a2, batch_res_DI_super, ...
    'VariableNames', {'Frequency_Hz', 'Lambda_m', ...
    'Source_X_m', 'Source_Y_m', 'Grid_Resolution_m', ...
    'Radial_Error_A1_m', 'Angular_Error_A1_deg', ...
    'Radial_Error_A2_m', 'Angular_Error_A2_deg', ...
    'Radial_Error_Super_m', 'Angular_Error_Super_deg', ...
    'Beamwidth_3dB_A1_deg', 'Beamwidth_3dB_A2_deg', 'Beamwidth_3dB_Super_deg', ...
    'Sidelobe_A1_dB', 'Sidelobe_A2_dB', 'Sidelobe_Super_dB', ...
    'DI_A1_dB', 'DI_A2_dB', 'DI_Super_dB'});

fprintf('\nResults table preview:\n');
disp(results_table);

csv_filename_out = fullfile(results_folder, 'freq_sweep_results.csv');
writetable(results_table, csv_filename_out);
fprintf('\nResults saved to: %s\n', csv_filename_out);

mat_filename_out = fullfile(results_folder, 'freq_sweep_results.mat');
save(mat_filename_out, 'results_table', 'batch_test_frequencies', ...
    'source_distance_m', 'source_angle_deg', 'delta_fixed', 'c_0', 'rho_0');
fprintf('Results also saved to: %s\n', mat_filename_out);


%% PUBLICATION FIGURES — FREQUENCY SWEEP RESULTS %%

col_a1 = [0.000, 0.447, 0.741]; %blue
col_a2 = [0.850, 0.325, 0.098]; %orange
col_super = [0.000, 0.000, 0.000]; %black
lw = 1.8;
ms = 7;
fn = 'Times New Roman';
fs_ax = 12;
fs_lab = 14;
fs_leg = 11;
freq_ticks = batch_test_frequencies;


% FIGURE: RADIAL ERROR VS FREQUENCY %

figure('Color', 'w', 'Position', [100 100 520 360]);
semilogy(batch_res_frequency_hz, batch_res_radial_error_a1, '-o', ...
    'Color', col_a1, 'LineWidth', lw, 'MarkerSize', ms, ...
    'MarkerFaceColor', col_a1, 'DisplayName', 'Array 1');
hold on;
semilogy(batch_res_frequency_hz, batch_res_radial_error_a2, '-s', ...
    'Color', col_a2, 'LineWidth', lw, 'MarkerSize', ms, ...
    'MarkerFaceColor', col_a2, 'DisplayName', 'Array 2');
semilogy(batch_res_frequency_hz, batch_res_radial_error_super, '-^', ...
    'Color', col_super, 'LineWidth', lw + 0.4, 'MarkerSize', ms, ...
    'MarkerFaceColor', col_super, 'DisplayName', 'Superimposed');
hold off;
xlabel('Frequency (Hz)', 'FontName', fn, 'FontSize', fs_lab);
ylabel('Radial error (m)', 'FontName', fn, 'FontSize', fs_lab);
legend('Location', 'best', 'FontName', fn, 'FontSize', fs_leg, 'Box', 'on');
set(gca, 'FontName', fn, 'FontSize', fs_ax, ...
    'XTick', freq_ticks, 'XTickLabelRotation', 0, ...
    'LineWidth', 1.0, 'Box', 'on', 'YGrid', 'on', 'XGrid', 'off');
exportgraphics(gcf, fullfile(results_folder, 'radial_error_vs_frequency.png'), ...
    'ContentType', 'image');


% FIGURE: ANGULAR ERROR VS FREQUENCY %

figure('Color', 'w', 'Position', [100 100 520 360]);
plot(batch_res_frequency_hz, batch_res_angular_error_a1, '-o', ...
    'Color', col_a1, 'LineWidth', lw, 'MarkerSize', ms, ...
    'MarkerFaceColor', col_a1, 'DisplayName', 'Array 1');
hold on;
plot(batch_res_frequency_hz, batch_res_angular_error_a2, '-s', ...
    'Color', col_a2, 'LineWidth', lw, 'MarkerSize', ms, ...
    'MarkerFaceColor', col_a2, 'DisplayName', 'Array 2');
plot(batch_res_frequency_hz, batch_res_angular_error_super, '-^', ...
    'Color', col_super, 'LineWidth', lw + 0.4, 'MarkerSize', ms, ...
    'MarkerFaceColor', col_super, 'DisplayName', 'Superimposed');
hold off;
xlabel('Frequency (Hz)', 'FontName', fn, 'FontSize', fs_lab);
ylabel('Angular error (deg)', 'FontName', fn, 'FontSize', fs_lab);
legend('Location', 'best', 'FontName', fn, 'FontSize', fs_leg, 'Box', 'on');
set(gca, 'FontName', fn, 'FontSize', fs_ax, ...
    'XTick', freq_ticks, 'LineWidth', 1.0, 'Box', 'on', ...
    'YGrid', 'on', 'XGrid', 'off');
exportgraphics(gcf, fullfile(results_folder, 'angular_error_vs_frequency.png'), ...
    'ContentType', 'image');


% FIGURE: BEAMWIDTH AND DIRECTIVITY INDEX (2x1 PANEL) %

figure('Color', 'w', 'Position', [100 100 520 600]);

subplot(2,1,1);
plot(batch_res_frequency_hz, batch_res_beamwidth_a1, '-o', ...
    'Color', col_a1, 'LineWidth', lw, 'MarkerSize', ms, ...
    'MarkerFaceColor', col_a1, 'DisplayName', 'Array 1');
hold on;
plot(batch_res_frequency_hz, batch_res_beamwidth_a2, '-s', ...
    'Color', col_a2, 'LineWidth', lw, 'MarkerSize', ms, ...
    'MarkerFaceColor', col_a2, 'DisplayName', 'Array 2');
plot(batch_res_frequency_hz, batch_res_beamwidth_super, '-^', ...
    'Color', col_super, 'LineWidth', lw + 0.4, 'MarkerSize', ms, ...
    'MarkerFaceColor', col_super, 'DisplayName', 'Superimposed');
hold off;
ylabel('$-3\,\mathrm{dB}$ beamwidth (deg)', 'Interpreter', 'latex', 'FontSize', fs_lab);
legend('Location', 'best', 'FontName', fn, 'FontSize', fs_leg, 'Box', 'on');
set(gca, 'FontName', fn, 'FontSize', fs_ax, ...
    'XTick', freq_ticks, 'XTickLabel', {}, ...
    'LineWidth', 1.0, 'Box', 'on', 'YGrid', 'on', 'XGrid', 'off');

subplot(2,1,2);
plot(batch_res_frequency_hz, batch_res_DI_a1, '-o', ...
    'Color', col_a1, 'LineWidth', lw, 'MarkerSize', ms, ...
    'MarkerFaceColor', col_a1, 'DisplayName', 'Array 1');
hold on;
plot(batch_res_frequency_hz, batch_res_DI_a2, '-s', ...
    'Color', col_a2, 'LineWidth', lw, 'MarkerSize', ms, ...
    'MarkerFaceColor', col_a2, 'DisplayName', 'Array 2');
plot(batch_res_frequency_hz, batch_res_DI_super, '-^', ...
    'Color', col_super, 'LineWidth', lw + 0.4, 'MarkerSize', ms, ...
    'MarkerFaceColor', col_super, 'DisplayName', 'Superimposed');
hold off;
xlabel('Frequency (Hz)', 'FontName', fn, 'FontSize', fs_lab);
ylabel('Directivity index (dB)', 'FontName', fn, 'FontSize', fs_lab);
legend('Location', 'best', 'FontName', fn, 'FontSize', fs_leg, 'Box', 'on');
set(gca, 'FontName', fn, 'FontSize', fs_ax, ...
    'XTick', freq_ticks, 'LineWidth', 1.0, 'Box', 'on', ...
    'YGrid', 'on', 'XGrid', 'off');
exportgraphics(gcf, fullfile(results_folder, 'beamwidth_DI_vs_frequency.png'), ...
    'ContentType', 'image');


%% FINAL SUMMARY %%

fprintf('\n========================================\n');
fprintf('<strong>FREQUENCY SWEEP TEST COMPLETE</strong>\n');
fprintf('========================================\n');

fprintf('\n%-12s %-16s %-16s %-16s\n', 'Freq (Hz)', 'A1 Radial (m)', ...
    'A2 Radial (m)', 'Super Radial (m)');
fprintf('%s\n', repmat('-', 1, 62));
for i = 1:batch_num_tests
    fprintf('%-12.0f %-16.4f %-16.4f %-16.4f\n', ...
        batch_res_frequency_hz(i), ...
        batch_res_radial_error_a1(i), ...
        batch_res_radial_error_a2(i), ...
        batch_res_radial_error_super(i));
end
fprintf('%s\n', repmat('-', 1, 62));

fprintf('\n%-12s %-16s %-16s %-16s\n', 'Freq (Hz)', 'A1 BW (deg)', ...
    'A2 BW (deg)', 'Super BW (deg)');
fprintf('%s\n', repmat('-', 1, 62));
for i = 1:batch_num_tests
    fprintf('%-12.0f %-16.1f %-16.1f %-16.1f\n', ...
        batch_res_frequency_hz(i), ...
        batch_res_beamwidth_a1(i), ...
        batch_res_beamwidth_a2(i), ...
        batch_res_beamwidth_super(i));
end
fprintf('%s\n', repmat('-', 1, 62));
fprintf('\nAll results saved to: %s\n\n', results_folder);


%% FUNCTION DEFINITIONS %%
% Replicated from the single-test AVA MVDR script to keep this script
% self-contained. Update both scripts if any function changes.


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
% Convert snapshots to frequency domain and compute cross-spectral matrix
function r_vs = create_vs_csm(snapshots_vs, bin_index, delta, bin_freqs, rho_0, N_v, num_bins, c_0)

    signal_dim = 3 * N_v;
    fx = fft(snapshots_vs, [], 2);
    fx = fx(:, bin_index, :);
    r_vs = zeros(signal_dim, signal_dim, num_bins);
    num_snap = size(snapshots_vs, 3);

    for jf = 1:num_bins
        freq = bin_freqs(jf);
        fx_bin = squeeze(fx(:, jf, :));
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

    plot(source_positions(:,1), source_positions(:,2), 'kx', ...
        'MarkerSize', 5, 'LineWidth', 1.5, 'DisplayName', 'True source');
    plot(all_est_positions(:,1), all_est_positions(:,2), 'ko', ...
        'MarkerSize', 10, 'LineWidth', 1.5, 'MarkerFaceColor', 'none', ...
        'DisplayName', 'Estimated source');

    ac = mean(vs_centres_arrays{array_num}(1:2,:), 2);
    plot(ac(1), ac(2), 'ks', 'MarkerSize', 9, 'LineWidth', 1, ...
        'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'DisplayName', 'Array');

    legend('Location', 'northeast', 'FontName', 'Times New Roman', ...
        'FontSize', 10, 'Color', 'w', 'TextColor', 'k', 'Box', 'on');
end


% FUNCTION: SUPERIMPOSED_LOCALISATION
% Mean-normalised product combination with ray intersection as source estimate
function [all_est_positions, max_db, response_db_norm_combined_out] = ...
    superimposed_localisation(response_db_combined, X_grid, Y_grid, y_scan, x_scan, ...
    num_sources, source_positions, N_a, grid_res, response_db_arrays, vs_centres_arrays)

    min_separation = max([0.25, 3 * grid_res]);

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

    % Ray intersection estimate from per-array peak bearings
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
                    fprintf('  Ray intersection outside grid — falling back to global peak\n');
                end
            else
                fprintf('  WARNING: array rays are parallel — falling back to global peak\n');
            end
        end
    end

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

    plot(source_positions(:,1), source_positions(:,2), 'kx', ...
        'MarkerSize', 5, 'LineWidth', 1.5, 'DisplayName', 'True source');
    plot(all_est_positions(:,1), all_est_positions(:,2), 'ko', ...
        'MarkerSize', 10, 'LineWidth', 1.5, 'MarkerFaceColor', 'none', ...
        'DisplayName', 'Estimated source');

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
% Near-field beam pattern for a single array
% Returns beam_metrics struct and raw per-angle linear response (mvdr_resp_linear)
% for use in compute_superimposed_beam_pattern
function [theta_deg, beam_pattern, fig_handle, beam_metrics, mvdr_resp_linear] = ...
    compute_nearfield_beam_pattern(r_vs, vs_centres, array_centre, source_positions, ...
    bin_freqs, c_0, rho_0, loading, num_bins, N_v, radius, plot_title)

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

    mvdr_resp_linear = sum(mvdr_responses, 2);
    beam_pattern = 10 * log10(mvdr_resp_linear + eps);
    beam_pattern = beam_pattern - max(beam_pattern);

    [fig_handle, beam_metrics] = plot_beam_pattern_polar(theta_rad, theta_deg, ...
        beam_pattern, source_positions, array_centre, radius, plot_title);

    print_beam_metrics(theta_deg, beam_pattern, beam_metrics);
end


% FUNCTION: COMPUTE_SUPERIMPOSED_BEAM_PATTERN
% Multiplies mean-normalised per-angle linear responses from all arrays
% before normalisation — simulation only, not valid for sequential recordings
% function [theta_deg, beam_pattern, fig_handle, beam_metrics] = ...
%     compute_superimposed_beam_pattern(mvdr_responses_per_array, system_centre, ...
%     source_positions, radius, plot_title)
% 
%     fprintf('\n<strong>Computing Superimposed Beam Pattern (simulation only):</strong>\n');
%     fprintf('  System centre: (%.3f, %.3f) m\n', system_centre);
%     fprintf('  Range: %.3f m\n', radius);
% 
%     N_arrays = length(mvdr_responses_per_array);
%     num_angles = length(mvdr_responses_per_array{1});
% 
%     theta_deg = linspace(0, 360, num_angles + 1);
%     theta_deg = theta_deg(1:end-1);
%     theta_rad = deg2rad(theta_deg);
% 
%     % Mean-normalise each array's response before multiplication to avoid
%     % absolute power differences between arrays biasing the result
%     super_linear = ones(num_angles, 1);
%     for array = 1:N_arrays
%         resp_linear = mvdr_responses_per_array{array};
%         resp_linear_norm = resp_linear / mean(resp_linear);
%         super_linear = super_linear .* resp_linear_norm;
%     end
% 
%     beam_pattern = 10 * log10(super_linear + eps);
%     beam_pattern = beam_pattern - max(beam_pattern);
% 
%     [fig_handle, beam_metrics] = plot_beam_pattern_polar(theta_rad, theta_deg, ...
%         beam_pattern, source_positions, system_centre, radius, plot_title);
% 
%     fprintf('  Array responses multiplied: %d arrays\n', N_arrays);
%     print_beam_metrics(theta_deg, beam_pattern, beam_metrics);
% end


% FUNCTION: PLOT_BEAM_PATTERN_POLAR
% Polar beam pattern plot — shared by individual and superimposed functions
% Returns fig_handle and beam_metrics struct
function [fig_handle, beam_metrics] = plot_beam_pattern_polar(theta_rad, theta_deg, ...
    beam_pattern, source_positions, array_centre, radius, plot_title)

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

    % Beam pattern metrics
    beam_metrics = struct();
    [~, peak_idx] = max(beam_pattern);

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
    end

    beam_metrics.sidelobe_level = NaN;
    bp_col = beam_pattern(:);
    if any(~above)
        beam_metrics.sidelobe_level = max(bp_col(~above));
    end

    beam_linear = 10.^(beam_pattern(:) / 10);
    integral_val = trapz(theta_rad(:), beam_linear);
    beam_metrics.directivity_index = 10 * log10(2 * pi / integral_val);
end


% FUNCTION: PRINT_BEAM_METRICS
% Print beam pattern statistics to command window
function print_beam_metrics(theta_deg, beam_pattern, beam_metrics)

    fprintf('\n<strong>Beam Pattern Statistics:</strong>\n');
    [~, peak_idx] = max(beam_pattern);
    fprintf('  Peak response: 0.00 dB at %.1f degrees\n', theta_deg(peak_idx));

    if ~isnan(beam_metrics.beamwidth_3dB)
        fprintf('  -3 dB beamwidth: %.1f degrees\n', beam_metrics.beamwidth_3dB);
    else
        fprintf('  WARNING: could not identify main lobe region\n');
    end
    if ~isnan(beam_metrics.sidelobe_level)
        fprintf('  Maximum sidelobe level: %.2f dB\n', beam_metrics.sidelobe_level);
    end
    fprintf('  Directivity index: %.2f dB\n', beam_metrics.directivity_index);
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


% FUNCTION: CREATE_VS_CSM_SIMULTANEOUS
% Full cross-spectral matrix for two AVS arrays (simulation only)
% Includes cross-array coherence terms — not valid for sequential recordings
function r_vs_full = create_vs_csm_simultaneous(snapshots_a1, snapshots_a2, ...
    bin_index, delta, bin_freqs, rho_0, N_v, num_bins, c_0)

    signal_dim_per_array = 3 * N_v;
    signal_dim_total = 2 * signal_dim_per_array;
    num_snap = size(snapshots_a1, 3);

    if size(snapshots_a2, 3) ~= num_snap
        error('Snapshot counts must match for simultaneous CSM construction.');
    end

    fx_a1 = fft(snapshots_a1, [], 2);
    fx_a2 = fft(snapshots_a2, [], 2);
    fx_a1 = fx_a1(:, bin_index, :);
    fx_a2 = fx_a2(:, bin_index, :);

    r_vs_full = zeros(signal_dim_total, signal_dim_total, num_bins);

    for jf = 1:num_bins
        freq = bin_freqs(jf);
        R_freq = zeros(signal_dim_total, signal_dim_total);

        for snap = 1:num_snap
            fx_snap_a1 = squeeze(fx_a1(:, jf, snap));
            fx_snap_a2 = squeeze(fx_a2(:, jf, snap));

            [p1, vx1, vy1] = vs_outputs_local(fx_snap_a1, delta, freq, rho_0, N_v, c_0);
            [p2, vx2, vy2] = vs_outputs_local(fx_snap_a2, delta, freq, rho_0, N_v, c_0);

            rho_c = rho_0 * c_0;
            vs_signal = zeros(signal_dim_total, 1);

            for n = 1:N_v
                idx1 = (n - 1) * 3 + (1:3);
                vs_signal(idx1) = [p1(n); rho_c * vx1(n); rho_c * vy1(n)];

                idx2 = signal_dim_per_array + (n - 1) * 3 + (1:3);
                vs_signal(idx2) = [p2(n); rho_c * vx2(n); rho_c * vy2(n)];
            end

            R_freq = R_freq + (vs_signal * vs_signal');
        end

        r_vs_full(:,:,jf) = R_freq / num_snap;
    end
end


% FUNCTION: VS_OUTPUTS_LOCAL
% Local copy of vs_outputs — keeps create_vs_csm_simultaneous self-contained
function [p_vs, vx_vs, vy_vs] = vs_outputs_local(tx_freq, delta, freq, rho_0, N_v, c_0)

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


% FUNCTION: VS_STEERING_VECTOR_COMBINED
% Augmented steering vector for two AVS arrays — concatenates per-array
% steering vectors to enable MVDR to exploit cross-array coherence in the full CSM
function v_combined = vs_steering_vector_combined(vs_centres_a1, vs_centres_a2, ...
    source_pos, freq, c_0, rho_0, N_v)

    v1 = vs_steering_vector_local(vs_centres_a1, source_pos, freq, c_0, rho_0, N_v);
    v2 = vs_steering_vector_local(vs_centres_a2, source_pos, freq, c_0, rho_0, N_v);

    v_combined = [v1; v2];
end


% FUNCTION: MVDR_VS_BEAMFORMING_SIMULTANEOUS
% MVDR beamforming using the full 6x6 CSM — exploits cross-array phase coherence
function response_db = mvdr_vs_beamforming_simultaneous(r_vs_full, ...
    vs_centres_a1, vs_centres_a2, candidate_points, bin_freqs, c_0, rho_0, ...
    loading, num_bins, N_v)

    num_points = size(candidate_points, 1);
    mvdr_responses = zeros(num_points, num_bins);

    for jf = 1:num_bins
        freq = bin_freqs(jf);
        rf = squeeze(r_vs_full(:,:,jf));

        [ur, er] = eig(rf);
        erv = real(diag(er));
        erv = max(erv, 1e-12);
        lambda = loading * max(erv);

        for n = 1:num_points
            source_pos = candidate_points(n,:).';
            v_combined = vs_steering_vector_combined(vs_centres_a1, vs_centres_a2, ...
                source_pos, freq, c_0, rho_0, N_v);

            rxv = (ur * diag(1 ./ (erv + lambda)) * ur') * v_combined;
            denominator = v_combined' * rxv;

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


% FUNCTION: VS_STEERING_VECTOR_LOCAL
% Local copy of vs_steering_vector — keeps simultaneous functions self-contained
function v_vs = vs_steering_vector_local(vs_centres, source_pos, freq, c_0, rho_0, N_v)

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
