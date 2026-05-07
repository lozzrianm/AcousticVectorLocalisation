clc; clear all; close all;

% Test C2 RESUME — picks up an interrupted C2 run without redoing
% the iterations that already have signal CSVs on disk.
%
% Usage:
%   1. Set `results_folder` below to the existing C2 output folder
%      (the one that already contains ula_*.csv and avs_*.csv files)
%   2. Run this script.
%
% What it does:
%   - For each (fi, ri, ai) combination in the sweep grid, checks
%     whether both ula_*.csv and avs_*.csv already exist in the
%     results folder.
%   - If both exist, loads them and recomputes the metrics. No new
%     signal generation. Fast.
%   - If either is missing, regenerates the missing signal(s) using
%     the same rng-seeded calls as the original script.
%   - At the end, writes c2_results.csv / .mat and produces all
%     crossover figures, identical to what the original script would
%     have produced if it had completed.
%
% Why this works:
%   - The original script writes per-run signal CSVs to disk during
%     the sweep, but only writes the aggregated results table and
%     figures at the very end. A mid-sweep timeout therefore leaves
%     all completed iteration's signal CSVs intact and recoverable.
%   - The MVDR processing, metric extraction, and figure generation
%     are deterministic given the same input signal CSV, so reloading
%     the cached signals and recomputing the metrics gives bit-exact
%     the same answers as a fresh run would have produced.
%
% Written by L Marshall


%% POINT THIS AT YOUR EXISTING C2 OUTPUT FOLDER %%

results_folder = 'crossarch_testC2_superimposedvsULA_2026-05-06_22-10-33';

if ~isfolder(results_folder)
    error('Results folder does not exist: %s\nEdit the results_folder line at the top of this script.', results_folder);
end


%% TEST PARAMETERS — must match BATCH_CROSSARCH_C2.m exactly %%

test_frequencies   = [630, 800, 1000, 1250, 1600, 2000, 2500];
distance_lambdas   = [0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0];
source_angles_deg  = [90, 110, 135, 160, 180];

num_freq  = length(test_frequencies);
num_dist  = length(distance_lambdas);
num_angle = length(source_angles_deg);
total_runs = num_freq * num_dist * num_angle;

c_0 = 340;
rho_0 = 1.02;

delta_avs = 0.04;
N_a = 2;
N_v = 1;
A1_centre = [0; 0; 0];
A2_centre = [0; -0.32; 0];

N_mics_ula = 8;
y_a_ula = -0.16;
total_aperture_ula = 0.36;
d_y_ula = total_aperture_ula / 7;

overlap = 0.5;
loading = 1e-4;
analysis_halfwidth_hz = 10;
target_binwidth_hz = 1;

grid_resolution_fixed = 0.005;
margin_fixed = 0.3;

ray_num_points = 200;
ray_range_lambdas = [0.05, 6.0];

rng_seed_base = 42;

fprintf('\n<strong>TEST C2 RESUME</strong>\n');
fprintf('Folder: %s\n', results_folder);
fprintf('Total runs to verify: %d\n\n', total_runs);


%% INITIALISE STORAGE — same shapes as the original script %%

res_freq          = zeros(total_runs, 1);
res_dist_lambda   = zeros(total_runs, 1);
res_angle_deg     = zeros(total_runs, 1);
res_lambda_m      = zeros(total_runs, 1);
res_source_x      = zeros(total_runs, 1);
res_source_y      = zeros(total_runs, 1);
res_radial_ula    = zeros(total_runs, 1);
res_radial_avs    = zeros(total_runs, 1);
res_angular_ula   = zeros(total_runs, 1);
res_angular_avs   = zeros(total_runs, 1);
res_bracket_ula   = zeros(total_runs, 1);
res_bracket_avs   = zeros(total_runs, 1);
res_sat_ula       = false(total_runs, 1);
res_sat_avs       = false(total_runs, 1);
res_beamwidth_ula = zeros(total_runs, 1);

radial_ula_tensor    = nan(num_freq, num_dist, num_angle);
radial_avs_tensor    = nan(num_freq, num_dist, num_angle);
angular_ula_tensor   = nan(num_freq, num_dist, num_angle);
angular_avs_tensor   = nan(num_freq, num_dist, num_angle);
bracket_ula_tensor   = nan(num_freq, num_dist, num_angle);
bracket_avs_tensor   = nan(num_freq, num_dist, num_angle);
sat_ula_tensor       = false(num_freq, num_dist, num_angle);
sat_avs_tensor       = false(num_freq, num_dist, num_angle);
beamwidth_ula_tensor = nan(num_freq, num_dist, num_angle);


%% SWEEP — reuse cached CSVs where present, regenerate where missing %%

run_idx = 0;
total_tic = tic;
n_cached = 0;
n_regen  = 0;

for fi = 1:num_freq
    test_freq = test_frequencies(fi);
    lambda    = c_0 / test_freq;

    for ri = 1:num_dist
        source_distance_m = distance_lambdas(ri) * lambda;

        for ai = 1:num_angle
            run_idx = run_idx + 1;
            theta_deg = source_angles_deg(ai);

            source_x = source_distance_m * cosd(theta_deg);
            source_y = source_distance_m * sind(theta_deg);


            % CHECK CACHE — do both signal CSVs already exist? %

            csv_ula = fullfile(results_folder, ...
                sprintf('ula_%04d_%dHz.csv', run_idx, round(test_freq)));
            csv_avs = fullfile(results_folder, ...
                sprintf('avs_%04d_%dHz.csv', run_idx, round(test_freq)));

            ula_exists = isfile(csv_ula);
            avs_exists = isfile(csv_avs);

            if ula_exists && avs_exists
                % Use cached signals
                n_cached = n_cached + 1;
                if mod(run_idx, 50) == 0 || run_idx <= 5
                    fprintf('[%d/%d] CACHED: f=%.0f Hz r=%.2f theta=%.0f\n', ...
                        run_idx, total_runs, test_freq, distance_lambdas(ri), theta_deg);
                end
            else
                % Regenerate missing signals using the same rng seed as
                % the original script
                n_regen = n_regen + 1;
                fprintf('[%d/%d] REGEN:  f=%.0f Hz r=%.2f theta=%.0f (ula=%d avs=%d)\n', ...
                    run_idx, total_runs, test_freq, distance_lambdas(ri), theta_deg, ...
                    ula_exists, avs_exists);

                seed = rng_seed_base + run_idx * 1000;

                if ~ula_exists
                    rng(seed);
                    batch_test_freq     = test_freq;
                    batch_test_source_x = source_x;
                    batch_test_source_y = source_y;
                    batch_csv_name      = csv_ula;
                    run('IDEAL_ULA_OUTPUT_BATCH_2904.m');
                end

                if ~avs_exists
                    rng(seed);
                    batch_test_freq     = test_freq;
                    batch_test_source_x = source_x;
                    batch_test_source_y = source_y;
                    batch_test_delta    = delta_avs;
                    batch_test_N_a      = 2;
                    batch_csv_name      = csv_avs;
                    run('SUPER_VA_OUTPUT_BATCH.m');
                end
            end


            % LOAD AND PROCESS — ULA %

            data_ula = readmatrix(csv_ula);
            time_ula = data_ula(:, 1);
            tx_ula = data_ula(:, 2:end).';

            n_arr = linspace((N_mics_ula - 1)/2, -(N_mics_ula - 1)/2, N_mics_ula);
            mic_positions_ula = [zeros(1, N_mics_ula);
                                 y_a_ula + d_y_ula * n_arr;
                                 zeros(1, N_mics_ula)];
            array_centre_ula = [0; y_a_ula];

            d_t_ula = time_ula(2) - time_ula(1);
            F_s_ula = 1 / d_t_ula;
            [bin_idx_ula, bin_freqs_ula, sfft_ula] = ...
                pick_bins(test_freq, F_s_ula, length(time_ula), ...
                          analysis_halfwidth_hz, target_binwidth_hz);
            window_ula = hanning(sfft_ula)';
            snap_ula = make_snapshots(tx_ula, sfft_ula, overlap, window_ula);
            r_ula = create_csm_pressure(snap_ula, bin_idx_ula, N_mics_ula, length(bin_freqs_ula));


            % LOAD AND PROCESS — DUAL AVS %

            data_avs_table = readtable(csv_avs);
            data_avs = table2array(data_avs_table);
            time_avs = data_avs(:, 1);
            N_ma = N_v * 4;

            tx_avs_arrays = cell(N_a, 1);
            vs_centres_arrays = cell(N_a, 1);
            mic_offsets = delta_avs / 2 * [-1, 1, 1, -1; -1, -1, 1, 1; 0, 0, 0, 0];

            for arr = 1:N_a
                col_start = 1 + (arr - 1) * N_ma + 1;
                col_end   = col_start + N_ma - 1;
                tx_avs_arrays{arr} = data_avs(:, col_start:col_end).';
                if arr == 1
                    centre = A1_centre;
                else
                    centre = A2_centre;
                end
                vs_centres_arrays{arr} = centre;
            end
            system_centre_avs = (A1_centre(1:2) + A2_centre(1:2)) / 2;

            d_t_avs = time_avs(2) - time_avs(1);
            F_s_avs = 1 / d_t_avs;
            [bin_idx_avs, bin_freqs_avs, sfft_avs] = ...
                pick_bins(test_freq, F_s_avs, length(time_avs), ...
                          analysis_halfwidth_hz, target_binwidth_hz);
            window_avs = hanning(sfft_avs)';

            r_avs_arrays = cell(N_a, 1);
            for arr = 1:N_a
                snap_arr = make_snapshots(tx_avs_arrays{arr}, sfft_avs, overlap, window_avs);
                r_avs_arrays{arr} = create_csm_avs(snap_arr, bin_idx_avs, ...
                    delta_avs, bin_freqs_avs, rho_0, N_v, length(bin_freqs_avs), c_0);
            end


            % SEARCH GRID %

            system_centre_y = (array_centre_ula(2) + system_centre_avs(2)) / 2;
            [x_scan, y_scan, X_grid, Y_grid, candidate_points, grid_res] = ...
                build_search_grid(source_distance_m, system_centre_y, grid_resolution_fixed, margin_fixed);


            % MVDR — ULA %

            response_ula_db = mvdr_pressure(r_ula, mic_positions_ula, ...
                candidate_points, bin_freqs_ula, c_0, loading);
            grid_ula = reshape(response_ula_db, length(y_scan), length(x_scan));
            [est_x_ula, est_y_ula] = find_peak_for_ula(grid_ula, x_scan, y_scan, ...
                source_x, source_y, theta_deg);


            % MVDR — superimposed dual AVS %

            response_avs_arrays_db = cell(N_a, 1);
            for arr = 1:N_a
                response_avs_arrays_db{arr} = mvdr_avs(r_avs_arrays{arr}, ...
                    vs_centres_arrays{arr}, candidate_points, bin_freqs_avs, ...
                    c_0, rho_0, loading, length(bin_freqs_avs), N_v);
            end

            ray_origins = zeros(2, N_a);
            ray_dirs = zeros(2, N_a);
            ray_ok = false(1, N_a);
            for arr = 1:N_a
                ac = vs_centres_arrays{arr}(1:2);
                arr_grid = reshape(response_avs_arrays_db{arr} - max(response_avs_arrays_db{arr}), ...
                    length(y_scan), length(x_scan));
                [px, py, ~] = refine_peak_2d(arr_grid, x_scan, y_scan);
                d = [px; py] - ac;
                if norm(d) > 1e-6
                    ray_origins(:, arr) = ac;
                    ray_dirs(:, arr) = d / norm(d);
                    ray_ok(arr) = true;
                end
            end

            est_x_avs = NaN; est_y_avs = NaN;
            valid_idx = find(ray_ok);
            if length(valid_idx) >= 2
                p1 = ray_origins(:, valid_idx(1)); d1 = ray_dirs(:, valid_idx(1));
                p2 = ray_origins(:, valid_idx(2)); d2 = ray_dirs(:, valid_idx(2));
                A = [d1, -d2]; b = p2 - p1;
                if abs(det(A)) > 1e-6
                    t = A \ b;
                    inter = p1 + t(1) * d1;
                    if inter(1) >= min(x_scan) && inter(1) <= max(x_scan) && ...
                       inter(2) >= min(y_scan) && inter(2) <= max(y_scan)
                        est_x_avs = inter(1);
                        est_y_avs = inter(2);
                    end
                end
            end

            super_lin = ones(size(response_avs_arrays_db{1}));
            for arr = 1:N_a
                lin_k = 10.^(response_avs_arrays_db{arr} / 10);
                super_lin = super_lin .* (lin_k / mean(lin_k));
            end
            super_db = 10 * log10(super_lin + eps);
            super_db = super_db - max(super_db);
            super_grid = reshape(super_db, length(y_scan), length(x_scan));

            if isnan(est_x_avs)
                [est_x_avs, est_y_avs, ~] = refine_peak_2d(super_grid, x_scan, y_scan);
            end


            % METRICS — radial and angular %

            radial_err_ula = norm([source_x - est_x_ula, source_y - est_y_ula]);
            radial_err_avs = norm([source_x - est_x_avs, source_y - est_y_avs]);

            true_b_ula = atan2(source_y - array_centre_ula(2), ...
                               source_x - array_centre_ula(1));
            est_b_ula  = atan2(est_y_ula - array_centre_ula(2), ...
                               est_x_ula - array_centre_ula(1));
            ang_err_ula = rad2deg(abs(wrap_to_pi(true_b_ula - est_b_ula)));

            true_b_avs = atan2(source_y - system_centre_avs(2), ...
                               source_x - system_centre_avs(1));
            est_b_avs  = atan2(est_y_avs - system_centre_avs(2), ...
                               est_x_avs - system_centre_avs(1));
            ang_err_avs = rad2deg(abs(wrap_to_pi(true_b_avs - est_b_avs)));


            % METRICS — -3 dB range bracket %

            bearing_ula = rad2deg(atan2(source_y - array_centre_ula(2), ...
                                        source_x - array_centre_ula(1)));
            bearing_avs = rad2deg(atan2(source_y - system_centre_avs(2), ...
                                        source_x - system_centre_avs(1)));

            ray_resp_ula = ray_mvdr_pressure(r_ula, mic_positions_ula, ...
                array_centre_ula, bearing_ula, ray_range_lambdas, ...
                ray_num_points, bin_freqs_ula, c_0, loading, lambda);
            distances_arr = linspace(ray_range_lambdas(1) * lambda, ...
                                     ray_range_lambdas(2) * lambda, ray_num_points);
            [~, bracket_ula, ~, ~, sat_ula] = analyse_ray_peak(distances_arr, ray_resp_ula);

            ray_resp_arrays_log = zeros(ray_num_points, N_a);
            for arr = 1:N_a
                ray_resp_arrays_log(:, arr) = ray_mvdr_avs(r_avs_arrays{arr}, ...
                    vs_centres_arrays{arr}, system_centre_avs, bearing_avs, ...
                    ray_range_lambdas, ray_num_points, bin_freqs_avs, c_0, ...
                    rho_0, loading, length(bin_freqs_avs), N_v, lambda);
            end
            super_lin_ray = ones(ray_num_points, 1);
            for arr = 1:N_a
                lin_k = 10.^(ray_resp_arrays_log(:, arr) / 10);
                super_lin_ray = super_lin_ray .* (lin_k / mean(lin_k));
            end
            ray_resp_super = 10 * log10(super_lin_ray + eps);
            ray_resp_super = ray_resp_super - max(ray_resp_super);
            [~, bracket_avs, ~, ~, sat_avs] = analyse_ray_peak(distances_arr, ray_resp_super);


            % METRICS — beamwidth (ULA only) %

            beamwidth_ula = beamwidth_at_range(r_ula, mic_positions_ula, ...
                array_centre_ula, source_distance_m, bin_freqs_ula, ...
                c_0, loading, 'pressure', N_v, rho_0);


            % STORE %

            res_freq(run_idx)         = test_freq;
            res_dist_lambda(run_idx)  = distance_lambdas(ri);
            res_angle_deg(run_idx)    = theta_deg;
            res_lambda_m(run_idx)     = lambda;
            res_source_x(run_idx)     = source_x;
            res_source_y(run_idx)     = source_y;
            res_radial_ula(run_idx)   = radial_err_ula;
            res_radial_avs(run_idx)   = radial_err_avs;
            res_angular_ula(run_idx)  = ang_err_ula;
            res_angular_avs(run_idx)  = ang_err_avs;
            res_bracket_ula(run_idx)  = bracket_ula;
            res_bracket_avs(run_idx)  = bracket_avs;
            res_sat_ula(run_idx)      = sat_ula;
            res_sat_avs(run_idx)      = sat_avs;
            res_beamwidth_ula(run_idx)= beamwidth_ula;

            radial_ula_tensor(fi, ri, ai)    = radial_err_ula;
            radial_avs_tensor(fi, ri, ai)    = radial_err_avs;
            angular_ula_tensor(fi, ri, ai)   = ang_err_ula;
            angular_avs_tensor(fi, ri, ai)   = ang_err_avs;
            bracket_ula_tensor(fi, ri, ai)   = bracket_ula;
            bracket_avs_tensor(fi, ri, ai)   = bracket_avs;
            sat_ula_tensor(fi, ri, ai)       = sat_ula;
            sat_avs_tensor(fi, ri, ai)       = sat_avs;
            beamwidth_ula_tensor(fi, ri, ai) = beamwidth_ula;

            close all;
        end
    end

    % Periodic snapshot save — write after each frequency completes so
    % a crash in mid-resume still leaves a partial results table on disk
    save(fullfile(results_folder, 'c2_results_partial.mat'), ...
        'res_freq', 'res_dist_lambda', 'res_angle_deg', 'res_lambda_m', ...
        'res_source_x', 'res_source_y', ...
        'res_radial_ula', 'res_radial_avs', 'res_angular_ula', 'res_angular_avs', ...
        'res_bracket_ula', 'res_bracket_avs', 'res_sat_ula', 'res_sat_avs', ...
        'res_beamwidth_ula', ...
        'radial_ula_tensor', 'radial_avs_tensor', ...
        'angular_ula_tensor', 'angular_avs_tensor', ...
        'bracket_ula_tensor', 'bracket_avs_tensor', ...
        'sat_ula_tensor', 'sat_avs_tensor', 'beamwidth_ula_tensor', ...
        'test_frequencies', 'distance_lambdas', 'source_angles_deg', ...
        'run_idx', 'fi');
    fprintf('  -> partial results snapshot saved after f=%.0f Hz (run %d)\n', test_freq, run_idx);
end

fprintf('\n<strong>C2 RESUME COMPLETE in %.1f min</strong>\n', toc(total_tic) / 60);
fprintf('Cached iterations: %d / %d (%.1f%% time saved)\n', ...
    n_cached, total_runs, 100 * n_cached / total_runs);
fprintf('Regenerated iterations: %d / %d\n\n', n_regen, total_runs);


%% SAVE FINAL %%

results_table = table(res_freq, res_dist_lambda, res_angle_deg, res_lambda_m, ...
    res_source_x, res_source_y, ...
    res_radial_ula, res_radial_avs, res_angular_ula, res_angular_avs, ...
    res_bracket_ula, res_bracket_avs, res_sat_ula, res_sat_avs, ...
    res_beamwidth_ula, ...
    'VariableNames', {'Frequency_Hz', 'Distance_Lambda', 'Angle_deg', 'Lambda_m', ...
    'Source_X_m', 'Source_Y_m', ...
    'Radial_ULA_m', 'Radial_DualAVS_m', 'Angular_ULA_deg', 'Angular_DualAVS_deg', ...
    'Bracket_ULA_m', 'Bracket_DualAVS_m', 'Saturated_ULA', 'Saturated_DualAVS', ...
    'Beamwidth_ULA_deg'});

writetable(results_table, fullfile(results_folder, 'c2_results.csv'));
save(fullfile(results_folder, 'c2_results.mat'), 'results_table', ...
    'radial_ula_tensor', 'radial_avs_tensor', ...
    'angular_ula_tensor', 'angular_avs_tensor', ...
    'bracket_ula_tensor', 'bracket_avs_tensor', ...
    'sat_ula_tensor', 'sat_avs_tensor', 'beamwidth_ula_tensor', ...
    'test_frequencies', 'distance_lambdas', 'source_angles_deg', ...
    'delta_avs', 'd_y_ula', 'N_mics_ula', 'A1_centre', 'A2_centre', 'c_0');


%% CROSSOVER FIGURES — identical to original script %%

set(groot, 'defaulttextinterpreter', 'latex');
set(groot, 'defaultaxesticklabelinterpreter', 'latex');
set(groot, 'defaultlegendinterpreter', 'latex');

for fi = 1:num_freq
    f = test_frequencies(fi);
    lam = c_0 / f;

    diff_radial = squeeze(radial_avs_tensor(fi,:,:) - radial_ula_tensor(fi,:,:)) / lam;
    plot_crossover_map(distance_lambdas, source_angles_deg, diff_radial, ...
        sprintf('Radial error advantage ($\\lambda$), $f = %.0f$ Hz', f), ...
        'Dual AVS optimal', 'ULA optimal', ...
        fullfile(results_folder, sprintf('c2_radial_crossover_%dHz', round(f))));

    diff_angular = squeeze(angular_avs_tensor(fi,:,:) - angular_ula_tensor(fi,:,:));
    plot_crossover_map(distance_lambdas, source_angles_deg, diff_angular, ...
        sprintf('Angular error advantage (deg), $f = %.0f$ Hz', f), ...
        'Dual AVS optimal', 'ULA optimal', ...
        fullfile(results_folder, sprintf('c2_angular_crossover_%dHz', round(f))));

    bracket_diff = squeeze(bracket_avs_tensor(fi,:,:) - bracket_ula_tensor(fi,:,:)) / lam;
    both_sat = squeeze(sat_ula_tensor(fi,:,:)) & squeeze(sat_avs_tensor(fi,:,:));
    bracket_diff(both_sat) = NaN;
    plot_crossover_map(distance_lambdas, source_angles_deg, bracket_diff, ...
        sprintf('$-3\\,\\mathrm{dB}$ range bracket advantage ($\\lambda$), $f = %.0f$ Hz', f), ...
        'Dual AVS optimal', 'ULA optimal', ...
        fullfile(results_folder, sprintf('c2_bracket_crossover_%dHz', round(f))));
end

fprintf('\nAll figures saved to: %s\n\n', results_folder);


%% LOCAL FUNCTIONS — copied verbatim from BATCH_CROSSARCH_C2.m %%

function [bin_index, bin_freqs, size_fft] = pick_bins(test_freq, F_s, N, halfwidth, target_bw)
    cycles_per_window = round(test_freq / target_bw);
    size_fft = round(cycles_per_window * F_s / test_freq);
    if size_fft > N / 4
        size_fft = 2^floor(log2(N / 4));
    end
    fft_vec = F_s * (0:(size_fft - 1)) / size_fft;
    f_lo = test_freq - halfwidth;
    f_hi = test_freq + halfwidth;
    bin_mask = (fft_vec >= f_lo) & (fft_vec <= f_hi);
    bin_index = find(bin_mask);
    bin_freqs = fft_vec(bin_index);
    if isempty(bin_index)
        [~, nearest] = min(abs(fft_vec(1:floor(size_fft/2)) - test_freq));
        bin_index = nearest;
        bin_freqs = fft_vec(nearest);
    end
end


function snapshots = make_snapshots(tx, size_fft, overlap, window)
    step = round(overlap * size_fft);
    num_snap = floor((size(tx, 2) - size_fft) / step) + 1;
    start_idx = 1 + (0:(num_snap - 1)) * step;
    idx_matrix = start_idx + (0:size_fft - 1)';
    snapshots = tx(:, idx_matrix);
    snapshots = reshape(snapshots, size(tx, 1), size_fft, num_snap);
    snapshots = snapshots .* reshape(window, [1, size_fft, 1]);
end


function r = create_csm_pressure(snapshots, bin_index, Nr, num_bins)
    fx = fft(snapshots, [], 2);
    fx = fx(:, bin_index, :);
    r = zeros(Nr, Nr, num_bins);
    for jf = 1:num_bins
        fx1 = squeeze(fx(:, jf, :));
        r(:,:,jf) = fx1 * fx1';
    end
    r = r / size(snapshots, 3);
end


function response_db = mvdr_pressure(r, mic_positions, candidate_points, bin_freqs, c_0, loading)
    num_points = size(candidate_points, 1);
    num_bins = length(bin_freqs);
    mvdr_responses = zeros(num_points, num_bins);
    for jf = 1:num_bins
        k = 2 * pi * bin_freqs(jf) / c_0;
        rf = squeeze(r(:,:,jf));
        [ur, er] = eig(rf);
        erv = max(real(diag(er)), 1e-12);
        for n = 1:num_points
            l_pt = candidate_points(n, :).';
            r_lm = sqrt(sum((mic_positions - l_pt).^2, 1));
            v = exp(-1i * k * r_lm).';
            cbf_out = real(v' * rf * v);
            lam = loading * cbf_out;
            if lam == 0, lam = loading; end
            rxv = (ur * diag(1 ./ (erv + lam)) * ur') * v;
            denom = v' * rxv;
            if abs(denom) > 1e-12
                vmvdr = rxv / denom;
                mvdr_responses(n, jf) = abs(vmvdr' * rf * vmvdr);
            end
        end
    end
    responses_sum = sum(mvdr_responses, 2);
    response_db = 10 * log10(responses_sum + eps);
    response_db = response_db - max(response_db);
end


function r_vs = create_csm_avs(snapshots_vs, bin_index, delta, bin_freqs, rho_0, N_v, num_bins, c_0)
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
            [p_vs, vx_vs, vy_vs] = vs_outputs(fx_snap, delta, freq, rho_0, N_v);
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


function [p_vs, vx_vs, vy_vs] = vs_outputs(tx_freq, delta, freq, rho_0, N_v)
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
        coeffs = M \ tx_freq(idx);
        p_vs(vs) = coeffs(1);
        vx_vs(vs) = -coeffs(2) / (1i * omega * rho_0);
        vy_vs(vs) = -coeffs(3) / (1i * omega * rho_0);
    end
end


function v_vs = vs_steering_spherical(vs_centres, source_pos, freq, c_0, rho_0, N_v)
    k = 2 * pi * freq / c_0;
    omega = 2 * pi * freq;
    v_vs = zeros(3 * N_v, 1);
    for n = 1:N_v
        r_vec = vs_centres(:, n) - source_pos;
        R = norm(r_vec);
        r_hat = r_vec / R;
        p_steer = exp(-1i * k * R) / R;
        common = exp(-1i * k * R) / (1i * omega * rho_0 * R^2);
        rho_c = rho_0 * c_0;
        vx_steer = rho_c * common * (1 + 1i * k * R) * r_hat(1);
        vy_steer = rho_c * common * (1 + 1i * k * R) * r_hat(2);
        idx = (n - 1) * 3 + (1:3);
        v_vs(idx) = [p_steer; vx_steer; vy_steer];
    end
end


function response_db = mvdr_avs(r_vs, vs_centres, candidate_points, bin_freqs, c_0, rho_0, loading, num_bins, N_v)
    num_points = size(candidate_points, 1);
    mvdr_responses = zeros(num_points, num_bins);
    for jf = 1:num_bins
        freq = bin_freqs(jf);
        rf = squeeze(r_vs(:,:,jf));
        [ur, er] = eig(rf);
        erv = max(real(diag(er)), 1e-12);
        lam = loading * max(erv);
        for n = 1:num_points
            source_pos = candidate_points(n, :).';
            v = vs_steering_spherical(vs_centres, source_pos, freq, c_0, rho_0, N_v);
            rxv = (ur * diag(1 ./ (erv + lam)) * ur') * v;
            denom = v' * rxv;
            if abs(denom) > 1e-12
                vmvdr = rxv / denom;
                mvdr_responses(n, jf) = abs(vmvdr' * rf * vmvdr);
            end
        end
    end
    responses_sum = sum(mvdr_responses, 2);
    response_db = 10 * log10(responses_sum + eps);
    response_db = response_db - max(response_db);
end


function [est_x, est_y] = find_peak_for_ula(grid_response, x_scan, y_scan, source_x, source_y, theta_deg)
    masked = grid_response;
    if source_x > 0
        masked(:, x_scan < 0) = -Inf;
    elseif source_x < 0
        masked(:, x_scan > 0) = -Inf;
    end
    [est_x, est_y, ~] = refine_peak_2d(masked, x_scan, y_scan);
end


function [x_scan, y_scan, X_grid, Y_grid, candidate_points, grid_res] = ...
    build_search_grid(source_distance, centre_y, grid_resolution_fixed, margin_fixed)
    max_extent = source_distance + margin_fixed;
    nx = round(2 * max_extent / grid_resolution_fixed) + 1;
    ny = round(2 * max_extent / grid_resolution_fixed) + 1;
    x_scan = linspace(-max_extent, max_extent, nx);
    y_scan = linspace(centre_y - max_extent, centre_y + max_extent, ny);
    [X_grid, Y_grid] = meshgrid(x_scan, y_scan);
    candidate_points = [X_grid(:), Y_grid(:), zeros(numel(X_grid), 1)];
    grid_res = mean([x_scan(2) - x_scan(1), y_scan(2) - y_scan(1)]);
end


function [est_x, est_y, refined_db] = refine_peak_2d(grid_response, x_scan, y_scan)
    [~, max_idx] = max(grid_response(:));
    [iy_pk, ix_pk] = ind2sub(size(grid_response), max_idx);
    est_x_disc = x_scan(ix_pk);
    est_y_disc = y_scan(iy_pk);
    dx = x_scan(2) - x_scan(1);
    dy = y_scan(2) - y_scan(1);
    if ix_pk < 2 || ix_pk > length(x_scan) - 1 || ...
       iy_pk < 2 || iy_pk > length(y_scan) - 1
        est_x = est_x_disc; est_y = est_y_disc;
        refined_db = grid_response(iy_pk, ix_pk);
        return;
    end
    fx_m = grid_response(iy_pk, ix_pk - 1);
    fx_0 = grid_response(iy_pk, ix_pk);
    fx_p = grid_response(iy_pk, ix_pk + 1);
    denom_x = fx_m - 2*fx_0 + fx_p;
    if abs(denom_x) > eps
        delta_ix = (fx_m - fx_p) / (2 * denom_x);
    else
        delta_ix = 0;
    end
    fy_m = grid_response(iy_pk - 1, ix_pk);
    fy_0 = grid_response(iy_pk, ix_pk);
    fy_p = grid_response(iy_pk + 1, ix_pk);
    denom_y = fy_m - 2*fy_0 + fy_p;
    if abs(denom_y) > eps
        delta_iy = (fy_m - fy_p) / (2 * denom_y);
    else
        delta_iy = 0;
    end
    delta_ix = max(-0.5, min(0.5, delta_ix));
    delta_iy = max(-0.5, min(0.5, delta_iy));
    est_x = est_x_disc + delta_ix * dx;
    est_y = est_y_disc + delta_iy * dy;
    refined_db = fx_0 - (fx_m - fx_p)^2 / (8 * denom_x);
end


function ray_response = ray_mvdr_pressure(r, mic_positions, array_centre, bearing_deg, ...
    range_lambdas, num_points, bin_freqs, c_0, loading, lambda)
    r_min = range_lambdas(1) * lambda;
    r_max = range_lambdas(2) * lambda;
    distances = linspace(r_min, r_max, num_points);
    rx = array_centre(1) + distances * cosd(bearing_deg);
    ry = array_centre(2) + distances * sind(bearing_deg);
    ray_points = [rx(:), ry(:), zeros(num_points, 1)];
    ray_response = mvdr_pressure(r, mic_positions, ray_points, bin_freqs, c_0, loading);
end


function ray_response = ray_mvdr_avs(r_vs, vs_centres, array_centre, bearing_deg, ...
    range_lambdas, num_points, bin_freqs, c_0, rho_0, loading, num_bins, N_v, lambda)
    r_min = range_lambdas(1) * lambda;
    r_max = range_lambdas(2) * lambda;
    distances = linspace(r_min, r_max, num_points);
    rx = array_centre(1) + distances * cosd(bearing_deg);
    ry = array_centre(2) + distances * sind(bearing_deg);
    ray_points = [rx(:), ry(:), zeros(num_points, 1)];
    ray_response = mvdr_avs(r_vs, vs_centres, ray_points, bin_freqs, c_0, rho_0, loading, num_bins, N_v);
end


function [peak_r, range_3dB, r_lower, r_upper, saturated] = analyse_ray_peak(distances, response_db)
    [pk_val, pk_idx] = max(response_db);
    if pk_idx < 2 || pk_idx > length(distances) - 1
        peak_r = distances(pk_idx);
    else
        dr = distances(2) - distances(1);
        f_m = response_db(pk_idx - 1);
        f_0 = response_db(pk_idx);
        f_p = response_db(pk_idx + 1);
        denom = f_m - 2 * f_0 + f_p;
        if abs(denom) > eps
            d_i = (f_m - f_p) / (2 * denom);
        else
            d_i = 0;
        end
        d_i = max(-0.5, min(0.5, d_i));
        peak_r = distances(pk_idx) + d_i * dr;
    end
    threshold = pk_val - 3;
    above = response_db(:) >= threshold;
    saturated = above(1) && above(end);
    if saturated
        range_3dB = distances(end) - distances(1);
        r_lower = distances(1);
        r_upper = distances(end);
    else
        d_above = diff([0; above; 0]);
        starts = find(d_above == 1);
        ends_arr = find(d_above == -1) - 1;
        region = find(starts <= pk_idx & ends_arr >= pk_idx, 1);
        if isempty(region)
            range_3dB = NaN; r_lower = NaN; r_upper = NaN;
        else
            r_lower = distances(starts(region));
            r_upper = distances(ends_arr(region));
            range_3dB = r_upper - r_lower;
        end
    end
end


function bw = beamwidth_at_range(r_in, geom, array_centre, radius, bin_freqs, ...
    c_0, loading, mode, N_v, rho_0)
    num_angles = 360;
    theta = linspace(0, 2*pi, num_angles + 1); theta = theta(1:end-1);
    cand = [array_centre(1) + radius * cos(theta);
            array_centre(2) + radius * sin(theta);
            zeros(1, num_angles)].';
    if strcmp(mode, 'pressure')
        resp = mvdr_pressure(r_in, geom, cand, bin_freqs, c_0, loading);
    else
        resp = mvdr_avs(r_in, geom, cand, bin_freqs, c_0, rho_0, loading, length(bin_freqs), N_v);
    end
    [pk_val, pk_idx] = max(resp);
    above = resp(:) >= (pk_val - 3);
    d_above = diff([0; above; 0]);
    starts = find(d_above == 1);
    ends_arr = find(d_above == -1) - 1;
    region = find(starts <= pk_idx & ends_arr >= pk_idx, 1);
    if isempty(region)
        bw = NaN;
    else
        a1 = rad2deg(theta(starts(region)));
        a2 = rad2deg(theta(ends_arr(region)));
        bw = a2 - a1;
        if bw < 0, bw = bw + 360; end
    end
end


function w = wrap_to_pi(x)
    w = mod(x + pi, 2*pi) - pi;
end


function plot_crossover_map(distance_lambdas, source_angles_deg, diff_data, ...
    cb_label, pos_label, neg_label, save_basename)
    plot_data = -diff_data;
    fn = 'Times New Roman';
    fs_ax = 13; fs_lab = 15;
    linthresh = 0.05;
    symlog = @(x) sign(x) .* log10(1 + abs(x) ./ linthresh);
    display_data = symlog(plot_data);
    cmax_data = max(abs(plot_data(:)), [], 'omitnan');
    if isnan(cmax_data) || cmax_data == 0
        cmax_data = 1;
    end
    cmax_disp = symlog(cmax_data);

    x_edges = compute_cell_edges(distance_lambdas);
    y_edges = compute_cell_edges(source_angles_deg);

    figure('Color', 'w', 'Position', [100 100 760 540]);
    data_padded = [display_data, nan(size(display_data, 1), 1); ...
                   nan(1, size(display_data, 2) + 1)];
    pc = pcolor(x_edges, y_edges, data_padded');
    pc.EdgeColor = 'none';
    axis xy;
    colormap(gca, redblue_colormap(256));
    clim([-cmax_disp, cmax_disp]);

    cb_real_ticks = [-5, -1, -0.2, 0, 0.2, 1, 5];
    cb_real_ticks = cb_real_ticks(abs(cb_real_ticks) <= cmax_data * 1.01);
    cb = colorbar;
    cb.Ticks = symlog(cb_real_ticks);
    cb.TickLabels = arrayfun(@(v) strip_zeros(sprintf('%.2f', v)), ...
        cb_real_ticks, 'UniformOutput', false);
    cb.Label.String = cb_label;
    cb.Label.Interpreter = 'latex';
    cb.Label.FontSize = fs_lab;
    cb.TickLabelInterpreter = 'latex';

    ax = gca;
    ax_pos = ax.Position;
    ax.Position = [ax_pos(1), ax_pos(2) + 0.04, ax_pos(3), ax_pos(4) - 0.08];
    cb.Position = [cb.Position(1), cb.Position(2) + 0.04, ...
                   cb.Position(3), cb.Position(4) - 0.08];
    cb_pos = cb.Position;
    annotation('textbox', ...
        [cb_pos(1) - 0.02, cb_pos(2) + cb_pos(4) - 0.01, ...
         cb_pos(3) + 0.08, 0.035], ...
        'String', pos_label, ...
        'Interpreter', 'latex', 'FontName', fn, 'FontSize', 11, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'EdgeColor', 'none', 'FitBoxToText', 'off');
    annotation('textbox', ...
        [cb_pos(1) - 0.02, cb_pos(2) - 0.037, ...
         cb_pos(3) + 0.08, 0.035], ...
        'String', neg_label, ...
        'Interpreter', 'latex', 'FontName', fn, 'FontSize', 11, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
        'EdgeColor', 'none', 'FitBoxToText', 'off');

    xlabel('Source distance $r/\lambda$', 'FontSize', fs_lab);
    ylabel('Source bearing (deg)', 'FontSize', fs_lab);

    x_label_strs = arrayfun(@(v) strip_zeros(sprintf('%.2f', v)), ...
        distance_lambdas, 'UniformOutput', false);
    y_label_strs = arrayfun(@(v) sprintf('%.0f', v), ...
        source_angles_deg, 'UniformOutput', false);

    set(gca, 'FontName', fn, 'FontSize', fs_ax, 'LineWidth', 1.0, ...
        'Box', 'on', 'TickLabelInterpreter', 'latex', ...
        'XTick', distance_lambdas, 'XTickLabel', x_label_strs, ...
        'YTick', source_angles_deg, 'YTickLabel', y_label_strs, ...
        'XMinorTick', 'on', 'YMinorTick', 'on', ...
        'TickDir', 'out', 'TickLength', [0.012 0.012]);
    ax = gca;
    ax.XAxis.MinorTickValues = x_edges;
    ax.YAxis.MinorTickValues = y_edges;
    xlim([x_edges(1), x_edges(end)]);
    ylim([y_edges(1), y_edges(end)]);

    save_thesis_figure(gcf, save_basename);
    close(gcf);
end


function save_thesis_figure(fig_handle, basename)
    exportgraphics(fig_handle, [basename, '.png'], 'ContentType', 'image', 'Resolution', 300);
    exportgraphics(fig_handle, [basename, '.pdf'], 'ContentType', 'vector');
    savefig(fig_handle, [basename, '.fig']);
end


function cmap = redblue_colormap(N)
    blue   = [0.10, 0.10, 1.00];
    white  = [1.00, 1.00, 1.00];
    orange = [0.9, 0.4, 0];
    if mod(N, 2) == 0
        half = N / 2;
        bot = [linspace(blue(1),  white(1), half)', ...
               linspace(blue(2),  white(2), half)', ...
               linspace(blue(3),  white(3), half)'];
        top = [linspace(white(1), orange(1), half)', ...
               linspace(white(2), orange(2), half)', ...
               linspace(white(3), orange(3), half)'];
    else
        half = (N - 1) / 2;
        bot = [linspace(blue(1),  white(1), half)', ...
               linspace(blue(2),  white(2), half)', ...
               linspace(blue(3),  white(3), half)'];
        top = [linspace(white(1), orange(1), half)', ...
               linspace(white(2), orange(2), half)', ...
               linspace(white(3), orange(3), half)'];
        bot = [bot; white];
    end
    cmap = [bot; top];
end


function edges = compute_cell_edges(test_values)
    v = test_values(:)';
    mids = (v(1:end-1) + v(2:end)) / 2;
    left_edge  = v(1)   - (mids(1)   - v(1));
    right_edge = v(end) + (v(end)    - mids(end));
    edges = [left_edge, mids, right_edge];
end


function s = strip_zeros(s)
    s = regexprep(s, '(\.\d*?)0+$', '$1');
    s = regexprep(s, '\.$', '');
end
