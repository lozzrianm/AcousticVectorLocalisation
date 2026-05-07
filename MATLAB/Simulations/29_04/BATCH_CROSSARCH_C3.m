clc; clear all; close all;

% Test C3 — AVS delta/lambda robustness check vs ULA reference
%
% Targeted line-plot variant of C1 to address the question: is the AVS's
% advantage in C1 robust to the choice of delta, or is it a cherry-picked
% optimum?
%
% Setup:
%   - One frequency (1000 Hz) — the operational midpoint of your test band
%   - One bearing (135 deg) — off-axis, the honest comparison case from Test B
%   - Two r/lambda points — one near-field (where AVS should dominate) and
%     one mid-field (where C1 should show the crossover)
%   - delta/lambda swept across the same Test B operational range
%   - AVS curves plotted vs delta/lambda; ULA shown as a horizontal line
%     (it doesn't depend on delta)
%
% Output figures (two panels each, radial and angular):
%   c3_radial_vs_delta.png/pdf/fig
%   c3_angular_vs_delta.png/pdf/fig
%   c3_bracket_vs_delta.png/pdf/fig
%
% Each panel shows AVS error (or bracket) as a function of delta/lambda
% with the ULA reference as a horizontal line. Crossover points where the
% AVS curve crosses the ULA line define the operational delta/lambda
% range over which the AVS is genuinely competitive at that geometry.
%
% Operational boundaries (from Test B / paper §3.4) marked as vertical lines:
%   delta/lambda < 0.07   — velocity SNR floor
%   delta/lambda > 0.33   — linearity breakdown (V_3 / V_6 degrade)
%   delta/lambda > 0.5    — spatial Nyquist (V_4 grating lobes)
%
% Written by L Marshall 06/05/2026


%% TEST PARAMETERS %%

% Single test frequency
test_freq = 1000;
c_0 = 340;
rho_0 = 1.02;
lambda = c_0 / test_freq;

% Single bearing — off-axis
source_angle_deg = 135;

% Two r/lambda points — one near, one mid
distance_lambdas = [0.3, 1.5];
num_dist = length(distance_lambdas);

% delta/lambda sweep — same as Test B, brackets the operational boundaries
delta_lambdas = [0.03, 0.05, 0.07, 0.10, 0.14, 0.18, 0.22, 0.27, 0.33, ...
                 0.40, 0.45];
num_delta = length(delta_lambdas);

total_runs = num_dist * num_delta;

% AVS configuration
N_a = 1;
N_v = 1;

% ULA configuration (matches IDEAL_ULA_OUTPUT_BATCH_2904.m)
N_mics_ula = 8;
y_a_ula = -0.16;
total_aperture_ula = 0.36;
d_y_ula = total_aperture_ula / 7;

% Processing
overlap = 0.5;
loading = 1e-4;
analysis_halfwidth_hz = 10;
target_binwidth_hz = 1;

grid_resolution_fixed = 0.005; %fixed grid step (m) — frequency-independent
margin_fixed = 0.3;

ray_num_points = 200;
ray_range_lambdas = [0.05, 6.0];

% Operational boundaries (from §3.4)
velocity_lower_bound = 0.07;
velocity_upper_bound = 0.33;
alias_threshold = 0.5;

rng_seed_base = 42;

timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
results_folder = sprintf('crossarch_C3_%s', timestamp);
mkdir(results_folder);

fprintf('\n<strong>TEST C3 — AVS delta/lambda ROBUSTNESS CHECK vs ULA</strong>\n');
fprintf('Frequency: %.0f Hz | Bearing: %.0f deg\n', test_freq, source_angle_deg);
fprintf('r/lambda: '); fprintf('%.2f ', distance_lambdas); fprintf('\n');
fprintf('delta/lambda (%d): ', num_delta); fprintf('%.2f  ', delta_lambdas); fprintf('\n');
fprintf('Total AVS runs: %d (plus %d ULA reference runs)\n', total_runs, num_dist);
fprintf('Results folder: %s\n\n', results_folder);


%% INITIALISE STORAGE %%

% AVS results: [num_dist x num_delta]
radial_avs    = nan(num_dist, num_delta);
angular_avs   = nan(num_dist, num_delta);
bracket_avs   = nan(num_dist, num_delta);
sat_avs       = false(num_dist, num_delta);
beamwidth_avs = nan(num_dist, num_delta);

% ULA reference: [num_dist x 1] (no delta dependence)
radial_ula    = nan(num_dist, 1);
angular_ula   = nan(num_dist, 1);
bracket_ula   = nan(num_dist, 1);
sat_ula       = false(num_dist, 1);
beamwidth_ula = nan(num_dist, 1);


%% RUN ULA REFERENCE (one per r/lambda) %%

fprintf('\n--- ULA reference runs ---\n');

for ri = 1:num_dist
    source_distance_m = distance_lambdas(ri) * lambda;
    source_x = source_distance_m * cosd(source_angle_deg);
    source_y = source_distance_m * sind(source_angle_deg);

    fprintf('\nULA: r/lambda=%.2f, source=(%.3f, %.3f) m\n', ...
        distance_lambdas(ri), source_x, source_y);

    seed = rng_seed_base + ri * 1000;
    rng(seed);
    batch_test_freq     = test_freq;
    batch_test_source_x = source_x;
    batch_test_source_y = source_y;
    batch_csv_name      = fullfile(results_folder, sprintf('ula_ref_r%.2f.csv', distance_lambdas(ri)));
    run('IDEAL_ULA_OUTPUT_BATCH_2904.m');

    data_ula = readmatrix(csv_filename);
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

    [x_scan, y_scan, ~, ~, candidate_points, ~] = ...
        build_search_grid(source_distance_m, ...
        array_centre_ula(2), grid_resolution_fixed, margin_fixed);


    response_ula_db = mvdr_pressure(r_ula, mic_positions_ula, ...
        candidate_points, bin_freqs_ula, c_0, loading);
    grid_ula = reshape(response_ula_db, length(y_scan), length(x_scan));
    [est_x, est_y] = find_peak_for_ula(grid_ula, x_scan, y_scan, ...
        source_x, source_y, source_angle_deg);

    radial_ula(ri) = norm([source_x - est_x, source_y - est_y]);
    true_b = atan2(source_y - array_centre_ula(2), source_x - array_centre_ula(1));
    est_b  = atan2(est_y - array_centre_ula(2), est_x - array_centre_ula(1));
    angular_ula(ri) = rad2deg(abs(wrap_to_pi(true_b - est_b)));

    bearing = rad2deg(true_b);
    ray_resp = ray_mvdr_pressure(r_ula, mic_positions_ula, array_centre_ula, ...
        bearing, ray_range_lambdas, ray_num_points, bin_freqs_ula, c_0, loading, lambda);
    distances_arr = linspace(ray_range_lambdas(1)*lambda, ray_range_lambdas(2)*lambda, ray_num_points);
    [~, bracket_ula(ri), ~, ~, sat_ula(ri)] = analyse_ray_peak(distances_arr, ray_resp);

    beamwidth_ula(ri) = beamwidth_at_range(r_ula, mic_positions_ula, ...
        array_centre_ula, source_distance_m, bin_freqs_ula, c_0, loading, ...
        'pressure', N_v, rho_0);

    fprintf('  Radial=%.4f m | Angular=%.3f deg | Bracket=%.4f m (sat=%d) | BW=%.1f deg\n', ...
        radial_ula(ri), angular_ula(ri), bracket_ula(ri), sat_ula(ri), beamwidth_ula(ri));

    close all;
end


%% RUN AVS SWEEP (delta/lambda x r/lambda) %%

fprintf('\n--- AVS sweep over delta/lambda ---\n');

run_idx = 0;
total_tic = tic;

for ri = 1:num_dist
    source_distance_m = distance_lambdas(ri) * lambda;
    source_x = source_distance_m * cosd(source_angle_deg);
    source_y = source_distance_m * sind(source_angle_deg);

    for di = 1:num_delta
        run_idx = run_idx + 1;
        delta_test = delta_lambdas(di) * lambda;

        fprintf('\n[%d/%d] r/lam=%.2f | delta/lam=%.2f (delta=%.4f m)\n', ...
            run_idx, total_runs, distance_lambdas(ri), delta_lambdas(di), delta_test);

        seed = rng_seed_base + ri * 1000 + di;
        rng(seed);
        batch_test_freq     = test_freq;
        batch_test_source_x = source_x;
        batch_test_source_y = source_y;
        batch_test_delta    = delta_test;
        batch_test_N_a      = 1;
        batch_csv_name      = fullfile(results_folder, ...
            sprintf('avs_r%.2f_d%.2f.csv', distance_lambdas(ri), delta_lambdas(di)));
        run('SUPER_VA_OUTPUT_BATCH.m');
        csv_avs = csv_filename;

        data_avs_table = readtable(csv_avs);
        data_avs = table2array(data_avs_table);
        time_avs = data_avs(:, 1);
        N_ma = N_v * 4;
        tx_avs = data_avs(:, 2:(N_ma + 1)).';

        A1_coord = [0; 0; 0];
        mic_offsets = delta_test / 2 * [-1, 1, 1, -1; -1, -1, 1, 1; 0, 0, 0, 0];
        vs_centres_avs = A1_coord;
        mic_positions_avs = vs_centres_avs + mic_offsets;
        array_centre_avs = [0; 0];

        d_t_avs = time_avs(2) - time_avs(1);
        F_s_avs = 1 / d_t_avs;
        [bin_idx_avs, bin_freqs_avs, sfft_avs] = ...
            pick_bins(test_freq, F_s_avs, length(time_avs), ...
                      analysis_halfwidth_hz, target_binwidth_hz);
        window_avs = hanning(sfft_avs)';
        snap_avs = make_snapshots(tx_avs, sfft_avs, overlap, window_avs);
        r_avs = create_csm_avs(snap_avs, bin_idx_avs, delta_test, ...
            bin_freqs_avs, rho_0, N_v, length(bin_freqs_avs), c_0);

        [x_scan, y_scan, ~, ~, candidate_points, ~] = ...
            build_search_grid(source_distance_m, ...
            array_centre_avs(2), grid_resolution_fixed, margin_fixed);

        response_avs_db = mvdr_avs(r_avs, vs_centres_avs, ...
            candidate_points, bin_freqs_avs, c_0, rho_0, loading, ...
            length(bin_freqs_avs), N_v);
        grid_avs = reshape(response_avs_db, length(y_scan), length(x_scan));
        [est_x, est_y, ~] = refine_peak_2d(grid_avs, x_scan, y_scan);

        radial_avs(ri, di) = norm([source_x - est_x, source_y - est_y]);
        true_b = atan2(source_y - array_centre_avs(2), source_x - array_centre_avs(1));
        est_b  = atan2(est_y - array_centre_avs(2), est_x - array_centre_avs(1));
        angular_avs(ri, di) = rad2deg(abs(wrap_to_pi(true_b - est_b)));

        bearing = rad2deg(true_b);
        ray_resp = ray_mvdr_avs(r_avs, vs_centres_avs, array_centre_avs, ...
            bearing, ray_range_lambdas, ray_num_points, bin_freqs_avs, ...
            c_0, rho_0, loading, length(bin_freqs_avs), N_v, lambda);
        distances_arr = linspace(ray_range_lambdas(1)*lambda, ray_range_lambdas(2)*lambda, ray_num_points);
        [~, bracket_avs(ri, di), ~, ~, sat_avs(ri, di)] = ...
            analyse_ray_peak(distances_arr, ray_resp);

        beamwidth_avs(ri, di) = beamwidth_at_range(r_avs, vs_centres_avs, ...
            array_centre_avs, source_distance_m, bin_freqs_avs, c_0, loading, ...
            'avs', N_v, rho_0);

        fprintf('  AVS: Radial=%.4f m | Angular=%.3f deg | Bracket=%.4f m (sat=%d) | BW=%.1f deg\n', ...
            radial_avs(ri, di), angular_avs(ri, di), bracket_avs(ri, di), ...
            sat_avs(ri, di), beamwidth_avs(ri, di));

        elapsed = toc(total_tic);
        eta = (elapsed / run_idx) * (total_runs - run_idx);
        fprintf('  Elapsed: %.1f min | ETA: %.1f min\n', elapsed/60, eta/60);

        close all;
    end
end

fprintf('\n<strong>C3 SWEEP COMPLETE in %.1f min</strong>\n', toc(total_tic) / 60);


%% SAVE RESULTS %%

% Save tensors and reference values to .mat for replotting
save(fullfile(results_folder, 'c3_results.mat'), ...
    'radial_avs', 'radial_ula', 'angular_avs', 'angular_ula', ...
    'bracket_avs', 'bracket_ula', 'sat_avs', 'sat_ula', ...
    'beamwidth_avs', 'beamwidth_ula', ...
    'delta_lambdas', 'distance_lambdas', 'test_freq', 'lambda', ...
    'source_angle_deg', 'velocity_lower_bound', 'velocity_upper_bound', ...
    'alias_threshold');

% CSV: long-form one row per (r/lambda, delta/lambda)
n_rows = num_dist * num_delta;
csv_dist  = zeros(n_rows, 1);
csv_delta = zeros(n_rows, 1);
csv_r_avs = zeros(n_rows, 1);
csv_r_ula = zeros(n_rows, 1);
csv_a_avs = zeros(n_rows, 1);
csv_a_ula = zeros(n_rows, 1);
csv_b_avs = zeros(n_rows, 1);
csv_b_ula = zeros(n_rows, 1);
csv_sat_avs = false(n_rows, 1);
csv_sat_ula = false(n_rows, 1);
csv_bw_avs = zeros(n_rows, 1);
csv_bw_ula = zeros(n_rows, 1);
ix = 0;
for ri = 1:num_dist
    for di = 1:num_delta
        ix = ix + 1;
        csv_dist(ix)    = distance_lambdas(ri);
        csv_delta(ix)   = delta_lambdas(di);
        csv_r_avs(ix)   = radial_avs(ri, di);
        csv_r_ula(ix)   = radial_ula(ri);
        csv_a_avs(ix)   = angular_avs(ri, di);
        csv_a_ula(ix)   = angular_ula(ri);
        csv_b_avs(ix)   = bracket_avs(ri, di);
        csv_b_ula(ix)   = bracket_ula(ri);
        csv_sat_avs(ix) = sat_avs(ri, di);
        csv_sat_ula(ix) = sat_ula(ri);
        csv_bw_avs(ix)  = beamwidth_avs(ri, di);
        csv_bw_ula(ix)  = beamwidth_ula(ri);
    end
end

results_table = table(csv_dist, csv_delta, csv_r_avs, csv_r_ula, ...
    csv_a_avs, csv_a_ula, csv_b_avs, csv_b_ula, ...
    csv_sat_avs, csv_sat_ula, csv_bw_avs, csv_bw_ula, ...
    'VariableNames', {'Distance_Lambda', 'Delta_Lambda', ...
    'Radial_AVS_m', 'Radial_ULA_m', 'Angular_AVS_deg', 'Angular_ULA_deg', ...
    'Bracket_AVS_m', 'Bracket_ULA_m', 'Saturated_AVS', 'Saturated_ULA', ...
    'Beamwidth_AVS_deg', 'Beamwidth_ULA_deg'});
writetable(results_table, fullfile(results_folder, 'c3_results.csv'));


%% LINE PLOTS — one figure per metric, two panels (one per r/lambda) %%

set(groot, 'defaulttextinterpreter', 'latex');
set(groot, 'defaultaxesticklabelinterpreter', 'latex');
set(groot, 'defaultlegendinterpreter', 'latex');

fn = 'Times New Roman';
fs_ax = 12; fs_lab = 14; fs_leg = 11;
col_avs = [0.85, 0.45, 0.10];
col_ula = [0.10, 0.30, 0.75];

% ---- Radial error vs delta/lambda ----
plot_metric_vs_delta(delta_lambdas, distance_lambdas, ...
    radial_avs / lambda, radial_ula / lambda, ...
    'Radial error ($r_\mathrm{err} / \lambda$)', ...
    velocity_lower_bound, velocity_upper_bound, alias_threshold, ...
    col_avs, col_ula, fn, fs_ax, fs_lab, fs_leg, 'log', ...
    fullfile(results_folder, 'c3_radial_vs_delta'));

% ---- Angular error vs delta/lambda ----
plot_metric_vs_delta(delta_lambdas, distance_lambdas, ...
    angular_avs, angular_ula, ...
    'Angular error (deg)', ...
    velocity_lower_bound, velocity_upper_bound, alias_threshold, ...
    col_avs, col_ula, fn, fs_ax, fs_lab, fs_leg, 'log', ...
    fullfile(results_folder, 'c3_angular_vs_delta'));

% ---- Range bracket vs delta/lambda ----
plot_metric_vs_delta(delta_lambdas, distance_lambdas, ...
    bracket_avs / lambda, bracket_ula / lambda, ...
    '$-3\,\mathrm{dB}$ range bracket ($\Delta r / \lambda$)', ...
    velocity_lower_bound, velocity_upper_bound, alias_threshold, ...
    col_avs, col_ula, fn, fs_ax, fs_lab, fs_leg, 'linear', ...
    fullfile(results_folder, 'c3_bracket_vs_delta'));

% ---- Beamwidth vs delta/lambda ----
plot_metric_vs_delta(delta_lambdas, distance_lambdas, ...
    beamwidth_avs, beamwidth_ula, ...
    '$-3\,\mathrm{dB}$ beamwidth (deg)', ...
    velocity_lower_bound, velocity_upper_bound, alias_threshold, ...
    col_avs, col_ula, fn, fs_ax, fs_lab, fs_leg, 'linear', ...
    fullfile(results_folder, 'c3_beamwidth_vs_delta'));

fprintf('\nAll C3 figures saved to: %s\n\n', results_folder);


%% PLOT FUNCTION %%

function plot_metric_vs_delta(delta_lambdas, distance_lambdas, ...
    avs_data, ula_data, ylabel_str, ...
    v_low, v_hi, alias_thr, col_avs, col_ula, ...
    fn, fs_ax, fs_lab, fs_leg, yscale, save_basename)
    % avs_data: [num_dist x num_delta]
    % ula_data: [num_dist x 1]
    % One panel per r/lambda

    num_dist = length(distance_lambdas);

    figure('Color', 'w', 'Position', [100 100 540 280 + 280 * num_dist]);

    for ri = 1:num_dist
        subplot(num_dist, 1, ri);

        % AVS curve
        if strcmp(yscale, 'log')
            semilogy(delta_lambdas, avs_data(ri, :), '-o', ...
                'Color', col_avs, 'LineWidth', 1.8, 'MarkerSize', 7, ...
                'MarkerFaceColor', col_avs, 'DisplayName', 'AVS');
        else
            plot(delta_lambdas, avs_data(ri, :), '-o', ...
                'Color', col_avs, 'LineWidth', 1.8, 'MarkerSize', 7, ...
                'MarkerFaceColor', col_avs, 'DisplayName', 'AVS');
        end
        hold on;

        % ULA reference horizontal line
        yline(ula_data(ri), '--', 'Color', col_ula, 'LineWidth', 2, ...
            'DisplayName', 'ULA (8 mics, 0.36 m)');

        % Operational boundary verticals
        xline(v_low,     'k:',  'LineWidth', 1.2, ...
            'DisplayName', 'Velocity SNR floor');
        xline(v_hi,      'k-.', 'LineWidth', 1.2, ...
            'DisplayName', 'Linearity breakdown');
        xline(alias_thr, 'k--', 'LineWidth', 1.2, ...
            'DisplayName', 'Spatial Nyquist');

        hold off;

        title(sprintf('$r/\\lambda = %.2f$', distance_lambdas(ri)), ...
            'FontSize', fs_lab);

        if ri == num_dist
            xlabel('Array spacing $\delta / \lambda$', 'FontSize', fs_lab);
        end
        ylabel(ylabel_str, 'FontSize', fs_lab);

        if ri == 1
            legend('Location', 'best', 'FontSize', fs_leg, 'Box', 'on', ...
                'NumColumns', 2);
        end

        set(gca, 'FontName', fn, 'FontSize', fs_ax, ...
            'XTick', delta_lambdas, 'LineWidth', 1.0, 'Box', 'on', ...
            'YGrid', 'on', 'XGrid', 'off', 'TickLabelInterpreter', 'latex');
    end

    save_thesis_figure(gcf, save_basename);
    close(gcf);
end


%% LOCAL FUNCTIONS — same as C1/C2 %%


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


function save_thesis_figure(fig_handle, basename)
    exportgraphics(fig_handle, [basename, '.png'], 'ContentType', 'image', 'Resolution', 300);
    exportgraphics(fig_handle, [basename, '.pdf'], 'ContentType', 'vector');
    savefig(fig_handle, [basename, '.fig']);
end
