clc; clear all; close all;

% PHASE SENSITIVITY MONTE CARLO EXPERIMENT — ULA vs AVA
%
% Compares the localisation robustness of three MVDR pipelines under
% per-channel time-delay calibration error:
%
%   (1) ULA-8        — ideal 8-mic single-recording ULA, full 8x8 CSM
%   (2) ULA-pooled   — two 4-mic recordings, block-diagonal pooled CSM
%                      (matches the experimental MEMS pipeline structure)
%   (3) AVA          — two single vector sensors, each a 4-mic colocated
%                      cluster (delta = 0.04 m), [P, Vx, Vy] LS-extracted
%                      and superimposed into a 6x6 CSM
%
% PHASE ERROR INJECTION
%   Per-channel time delays drawn from N(0, sigma_tau^2). Applied as
%   exp(-1i*2*pi*f*tau_n) on frequency-domain snapshots before CSM
%   accumulation (and before vs_outputs LS-regression for AVA). sigma_tau
%   chosen so 2*pi*f*sigma_tau corresponds to a target RMS phase error in
%   degrees at the test frequency.
%
% MONTE CARLO
%   For each phase RMS level, n_trials independent draws of the 8-element
%   tau vector. The SAME tau vector is applied to all three pipelines
%   within a trial (sub-array assignment: tau(1:4) -> sub1, tau(5:8) -> sub2).
%   This isolates between-pipeline differences from between-realisation
%   differences.
%
% OUTPUT
%   Median and 25/75 percentile bands of radial localisation error vs
%   phase error level, per pipeline. CSV summary and PDF figure.
%
% RUNTIME
%   With defaults (50 trials x 8 levels, ~4900 grid points), expect
%   ~10-20 min on a desktop machine. Reduce n_trials or shrink the
%   search grid for faster iteration during development.
%
% Written by L Marshall, structure adapted from IN26 batch scripts.

%% EXPERIMENT PARAMETERS %%

% Test frequency (single tone)
f = 1000;             %Hz
omega = 2 * pi * f;
c_0 = 340;            %speed of sound (m/s)
rho_0 = 1.02;         %air density (kg/m^3)
lambda = c_0 / f;
k = 2 * pi / lambda;

% Source position
src_x = -0.4;
src_y = 0.0;
src_z = 0.0;
src = [src_x; src_y; src_z];

% Phase error sweep — RMS phase error (deg) at test frequency
phase_levels_deg = [0.0, 0.5, 1.0, 2.0, 5.0, 10.0, 12.0, 20.0];
n_levels = length(phase_levels_deg);

% Monte Carlo trials per level
n_trials = 50;

% Synthesis parameters
num_snap = 200;       %frequency-domain snapshots per trial
loading = 1e-4;       %MVDR adaptive regularisation
Q_o = 0.01;           %source amplitude

% Search grid — fixed across all levels and pipelines for fair comparison
grid_pts_per_lambda = 20;
grid_margin = 0.2;    %physical margin around source distance (m)
grid_res = lambda / grid_pts_per_lambda;
grid_extent = norm([src_x, src_y]) + grid_margin;
y_centre_grid = -0.16;  %geometric centre of the longest array

x_scan = -grid_extent:grid_res:grid_extent;
y_scan = (y_centre_grid - grid_extent):grid_res:(y_centre_grid + grid_extent);
[X_grid, Y_grid] = meshgrid(x_scan, y_scan);
grid_xyz = [X_grid(:), Y_grid(:), zeros(numel(X_grid), 1)].';
n_grid = size(grid_xyz, 2);

% Output folder
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
results_folder = sprintf('phase_sensitivity_%s', timestamp);
mkdir(results_folder);

fprintf('\n<strong>PHASE SENSITIVITY MONTE CARLO</strong>\n');
fprintf('Source: (%.3f, %.3f) m at %.0f Hz | lambda = %.4f m\n', src_x, src_y, f, lambda);
fprintf('Phase RMS levels (deg): '); fprintf('%.1f ', phase_levels_deg); fprintf('\n');
fprintf('Trials per level: %d | Snapshots per trial: %d\n', n_trials, num_snap);
fprintf('Grid: %d x %d points, resolution %.4f m\n', length(x_scan), length(y_scan), grid_res);
fprintf('Output folder: %s\n\n', results_folder);


%% ARRAY GEOMETRIES %%

% ULA-8 (single recording)
d_y_ula = 0.36 / 7;
y_a_ula8 = -0.16;
n_vec_8 = linspace(3.5, -3.5, 8);
ula8_pos = [zeros(1, 8); y_a_ula8 + d_y_ula * n_vec_8; zeros(1, 8)];

% ULA-pooled sub-arrays (4 mics each)
y_centre_a1 = 0.02 - 1.5 * d_y_ula;
y_centre_a2 = -0.34 + 1.5 * d_y_ula;
n_vec_4 = linspace(1.5, -1.5, 4);
ulap_a1_pos = [zeros(1, 4); y_centre_a1 + d_y_ula * n_vec_4; zeros(1, 4)];
ulap_a2_pos = [zeros(1, 4); y_centre_a2 + d_y_ula * n_vec_4; zeros(1, 4)];
ulap_pos_full = [ulap_a1_pos, ulap_a2_pos];

% AVA — two VS clusters at sub-array centres, delta = 0.04 m
delta = 0.04;
ava_centres = [0, 0; 0, -0.32; 0, 0];   %(3, 2)
mic_offsets = [-delta/2, +delta/2, +delta/2, -delta/2;
               -delta/2, -delta/2, +delta/2, +delta/2;
               0, 0, 0, 0];
ava_a1_pos = ava_centres(:, 1) + mic_offsets;
ava_a2_pos = ava_centres(:, 2) + mic_offsets;

fprintf('Array geometries:\n');
fprintf('  ULA-8       aperture %.4f m, centre y = %.4f m\n', ...
    max(ula8_pos(2,:)) - min(ula8_pos(2,:)), mean(ula8_pos(2,:)));
fprintf('  ULA-pooled  aperture %.4f m, sub-array centres y = %.4f, %.4f m\n', ...
    max([ulap_a1_pos(2,:) ulap_a2_pos(2,:)]) - min([ulap_a1_pos(2,:) ulap_a2_pos(2,:)]), ...
    y_centre_a1, y_centre_a2);
fprintf('  AVA         VS centres y = %.4f, %.4f m, delta = %.3f m\n\n', ...
    ava_centres(2, 1), ava_centres(2, 2), delta);


%% RESULTS STORAGE %%

err_ula8   = nan(n_levels, n_trials);
err_pooled = nan(n_levels, n_trials);
err_ava    = nan(n_levels, n_trials);


%% MAIN LOOP %%

% Set RNG so results are reproducible across runs of this script
rng(42);

for li = 1:n_levels
    sigma_deg = phase_levels_deg(li);
    sigma_phi_rad = deg2rad(sigma_deg);
    sigma_tau = sigma_phi_rad / (2 * pi * f);

    fprintf('Level %d/%d: %.1f deg RMS (sigma_tau = %.3f us)\n', ...
        li, n_levels, sigma_deg, sigma_tau * 1e6);

    for ti = 1:n_trials
        % Draw 8 per-channel time delays for this trial
        taus = sigma_tau * randn(8, 1);

        % Generate clean frequency-domain snapshots (independent realisations
        % for each pipeline so the snapshots are not artificially correlated
        % across pipelines — only the tau vector is shared)
        P_ula8  = synth_snapshots(ula8_pos,    src, num_snap, f, c_0, rho_0, Q_o);
        P_ulap1 = synth_snapshots(ulap_a1_pos, src, num_snap, f, c_0, rho_0, Q_o);
        P_ulap2 = synth_snapshots(ulap_a2_pos, src, num_snap, f, c_0, rho_0, Q_o);
        P_ava1  = synth_snapshots(ava_a1_pos,  src, num_snap, f, c_0, rho_0, Q_o);
        P_ava2  = synth_snapshots(ava_a2_pos,  src, num_snap, f, c_0, rho_0, Q_o);

        % Inject phase error — same tau vector across pipelines, split by sub-array
        P_ula8  = apply_phase_err(P_ula8,  taus,     f);
        P_ulap1 = apply_phase_err(P_ulap1, taus(1:4), f);
        P_ulap2 = apply_phase_err(P_ulap2, taus(5:8), f);
        P_ava1  = apply_phase_err(P_ava1,  taus(1:4), f);
        P_ava2  = apply_phase_err(P_ava2,  taus(5:8), f);

        % --- ULA-8 ---
        R8 = (P_ula8 * P_ula8') / num_snap;
        e8 = mvdr_localise_pressure(R8, ula8_pos, grid_xyz, k, loading);
        err_ula8(li, ti) = norm(e8(1:2) - src(1:2));

        % --- ULA-pooled ---
        Rp1 = (P_ulap1 * P_ulap1') / num_snap;
        Rp2 = (P_ulap2 * P_ulap2') / num_snap;
        Rp = zeros(8, 8);
        Rp(1:4, 1:4) = Rp1;
        Rp(5:8, 5:8) = Rp2;
        ep = mvdr_localise_pressure(Rp, ulap_pos_full, grid_xyz, k, loading);
        err_pooled(li, ti) = norm(ep(1:2) - src(1:2));

        % --- AVA ---
        vs1 = vs_outputs(P_ava1, delta, omega, rho_0, c_0);
        vs2 = vs_outputs(P_ava2, delta, omega, rho_0, c_0);
        sig_avs = [vs1; vs2];               %(6, num_snap)
        Rava = (sig_avs * sig_avs') / num_snap;
        ea = mvdr_localise_avs(Rava, ava_centres, grid_xyz, k, loading);
        err_ava(li, ti) = norm(ea(1:2) - src(1:2));
    end

    % Quick progress summary at this level
    fprintf('  ULA-8:      median %.4f m  IQR [%.4f, %.4f]\n', ...
        median(err_ula8(li, :)), prctile(err_ula8(li, :), 25), prctile(err_ula8(li, :), 75));
    fprintf('  ULA-pooled: median %.4f m  IQR [%.4f, %.4f]\n', ...
        median(err_pooled(li, :)), prctile(err_pooled(li, :), 25), prctile(err_pooled(li, :), 75));
    fprintf('  AVA:        median %.4f m  IQR [%.4f, %.4f]\n', ...
        median(err_ava(li, :)),    prctile(err_ava(li, :), 25),    prctile(err_ava(li, :), 75));
end


%% SUMMARY STATISTICS %%

med_8   = median(err_ula8, 2);
p25_8   = prctile(err_ula8, 25, 2);
p75_8   = prctile(err_ula8, 75, 2);

med_p   = median(err_pooled, 2);
p25_p   = prctile(err_pooled, 25, 2);
p75_p   = prctile(err_pooled, 75, 2);

med_a   = median(err_ava, 2);
p25_a   = prctile(err_ava, 25, 2);
p75_a   = prctile(err_ava, 75, 2);


%% RESULTS TABLE %%

results_table = table(phase_levels_deg(:), ...
    med_8, p25_8, p75_8, ...
    med_p, p25_p, p75_p, ...
    med_a, p25_a, p75_a, ...
    'VariableNames', {'PhaseRMS_deg', ...
        'ULA8_median_m', 'ULA8_p25_m', 'ULA8_p75_m', ...
        'Pooled_median_m', 'Pooled_p25_m', 'Pooled_p75_m', ...
        'AVA_median_m', 'AVA_p25_m', 'AVA_p75_m'});

fprintf('\n<strong>RESULTS TABLE</strong>\n');
disp(results_table);

csv_out = fullfile(results_folder, 'phase_sensitivity_results.csv');
writetable(results_table, csv_out);
fprintf('Saved CSV: %s\n', csv_out);

% Also save raw per-trial errors for replotting
mat_out = fullfile(results_folder, 'phase_sensitivity_raw.mat');
save(mat_out, 'phase_levels_deg', 'n_trials', 'err_ula8', 'err_pooled', 'err_ava', ...
    'f', 'src', 'num_snap', 'loading');
fprintf('Saved raw data: %s\n\n', mat_out);


%% MAIN FIGURE — PERCENTILE BANDS %%

col_ula8   = [0.000, 0.447, 0.741];   %blue
col_pooled = [0.466, 0.674, 0.188];   %green
col_ava    = [0.850, 0.325, 0.098];   %orange-red

figure('Color', 'w', 'Position', [100 100 720 480]);
hold on;

% Shaded percentile bands
fill_band(phase_levels_deg, p25_8, p75_8, col_ula8, 0.15);
fill_band(phase_levels_deg, p25_p, p75_p, col_pooled, 0.15);
fill_band(phase_levels_deg, p25_a, p75_a, col_ava, 0.15);

% Median lines
h_8 = plot(phase_levels_deg, med_8, '-o', 'Color', col_ula8, ...
    'LineWidth', 2.0, 'MarkerSize', 7, 'MarkerFaceColor', col_ula8);
h_p = plot(phase_levels_deg, med_p, '-s', 'Color', col_pooled, ...
    'LineWidth', 2.0, 'MarkerSize', 7, 'MarkerFaceColor', col_pooled);
h_a = plot(phase_levels_deg, med_a, '-^', 'Color', col_ava, ...
    'LineWidth', 2.0, 'MarkerSize', 7, 'MarkerFaceColor', col_ava);

% Reference: grid resolution floor
yline(grid_res, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.0, 'HandleVisibility', 'off');
text(phase_levels_deg(end), grid_res, '  grid res', ...
    'FontName', 'Times New Roman', 'FontSize', 9, ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Color', [0.4 0.4 0.4]);

set(gca, 'YScale', 'log', 'FontName', 'Times New Roman', 'FontSize', 12, ...
    'LineWidth', 1.0, 'Box', 'on', 'YGrid', 'on', 'XGrid', 'off');
xlabel('RMS phase error per channel at $f$ (deg)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Radial localisation error (m)', 'Interpreter', 'latex', 'FontSize', 14);
title(sprintf('Phase Sensitivity: $f = %d$ Hz, source at $(%.2f, %.2f)$ m', ...
    f, src_x, src_y), 'Interpreter', 'latex', 'FontSize', 13);
legend([h_8, h_p, h_a], {'ULA-8 (full CSM)', 'ULA-pooled (block-diag CSM)', 'AVA (2x VS)'}, ...
    'Location', 'best', 'FontName', 'Times New Roman', 'FontSize', 11, 'Box', 'on');

exportgraphics(gcf, fullfile(results_folder, 'phase_sensitivity.pdf'), 'ContentType', 'vector');
exportgraphics(gcf, fullfile(results_folder, 'phase_sensitivity.png'), 'Resolution', 220);
fprintf('Saved figure: phase_sensitivity.pdf / .png\n');


%% SECONDARY FIGURE — ALL TRIALS SCATTER %%
% Useful for seeing bimodal failure modes (the "either it works or it jumps
% to a sidelobe" behaviour that shifts only the upper percentile)

figure('Color', 'w', 'Position', [100 100 720 480]);
hold on;
for li = 1:n_levels
    x_jit = phase_levels_deg(li) + 0.05 * randn(1, n_trials);
    scatter(x_jit, err_ula8(li, :), 18, col_ula8, 'filled', 'MarkerFaceAlpha', 0.4);
    scatter(x_jit, err_pooled(li, :), 18, col_pooled, 'filled', 'MarkerFaceAlpha', 0.4);
    scatter(x_jit, err_ava(li, :), 18, col_ava, 'filled', 'MarkerFaceAlpha', 0.4);
end
plot(phase_levels_deg, med_8, '-', 'Color', col_ula8, 'LineWidth', 2);
plot(phase_levels_deg, med_p, '-', 'Color', col_pooled, 'LineWidth', 2);
plot(phase_levels_deg, med_a, '-', 'Color', col_ava, 'LineWidth', 2);

set(gca, 'YScale', 'log', 'FontName', 'Times New Roman', 'FontSize', 12, ...
    'LineWidth', 1.0, 'Box', 'on', 'YGrid', 'on');
xlabel('RMS phase error per channel at $f$ (deg)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Radial error per trial (m)', 'Interpreter', 'latex', 'FontSize', 14);
title('Per-trial scatter (jittered) — exposes bimodal failure', ...
    'Interpreter', 'latex', 'FontSize', 12);

% Manual legend via dummy scatter handles (transparency confuses default legend)
hold on;
h1 = scatter(NaN, NaN, 30, col_ula8, 'filled');
h2 = scatter(NaN, NaN, 30, col_pooled, 'filled');
h3 = scatter(NaN, NaN, 30, col_ava, 'filled');
legend([h1, h2, h3], {'ULA-8', 'ULA-pooled', 'AVA'}, ...
    'Location', 'best', 'FontName', 'Times New Roman', 'FontSize', 11);

exportgraphics(gcf, fullfile(results_folder, 'phase_sensitivity_scatter.pdf'), 'ContentType', 'vector');
exportgraphics(gcf, fullfile(results_folder, 'phase_sensitivity_scatter.png'), 'Resolution', 220);
fprintf('Saved figure: phase_sensitivity_scatter.pdf / .png\n\n');

fprintf('<strong>EXPERIMENT COMPLETE</strong>\n');
fprintf('All outputs in: %s\n\n', results_folder);


%% ============================ FUNCTIONS ============================ %%


% FUNCTION: SYNTH_SNAPSHOTS
% Generate frequency-domain pressure snapshots at a single frequency bin f
% for a monopole-derivative source. Returns (N_mic x num_snap) complex.
% Each snapshot has an independent random source phase so snapshots are
% decorrelated and the CSM has the expected rank-1 + noise structure.
function P = synth_snapshots(mic_xyz, src_xyz, num_snap, f, c_0, rho_0, Q_o)
    N_mic = size(mic_xyz, 2);
    omega = 2 * pi * f;
    k = omega / c_0;
    R = sqrt(sum((mic_xyz - src_xyz).^2, 1)).';   %(N_mic x 1)
    amp = (rho_0 ./ (4 * pi * R)) * omega;        %monopole derivative magnitude
    src_phase = 2 * pi * rand(1, num_snap);       %random source phase per snapshot
    src_complex = Q_o * exp(1i * src_phase);      %(1 x num_snap)
    P = (amp .* exp(-1i * k * R)) * src_complex;  %(N_mic x num_snap)
end


% FUNCTION: APPLY_PHASE_ERR
% Apply per-channel time delay tau as exp(-1i*2*pi*f*tau) to each row of
% frequency-domain snapshots. Matches the convention used in the
% experimental scripts (note: see calibration write-up re sign convention).
function P_err = apply_phase_err(P, taus, f)
    corr = exp(-1i * 2 * pi * f * taus(:));  %(N_mic x 1)
    P_err = P .* corr;                       %broadcast to (N_mic x num_snap)
end


% FUNCTION: VS_OUTPUTS
% Reproduce the AVS LS pressure-velocity decomposition from
% IN26_AVA_Experimental_SUPER_MVDR.m line 908. For each VS cluster (4 mics
% at +/-delta/2 corners), solves the 4x3 system M*[P;Vx;Vy] = p_vector for
% each snapshot. Returns [P, rho*c*Vx, rho*c*Vy] stacked as 3xnum_snap.
function vs = vs_outputs(P_freq, delta, omega, rho_0, c_0)
    mic_pos_2d = [-delta/2, -delta/2;
                  +delta/2, -delta/2;
                  +delta/2, +delta/2;
                  -delta/2, +delta/2];
    M = [ones(4, 1), mic_pos_2d];            %(4 x 3)
    coeffs = M \ P_freq;                     %(3 x num_snap) least-squares
    P_vs  = coeffs(1, :);
    Vx_vs = -coeffs(2, :) / (1i * omega * rho_0);
    Vy_vs = -coeffs(3, :) / (1i * omega * rho_0);
    rho_c = rho_0 * c_0;
    vs = [P_vs; rho_c * Vx_vs; rho_c * Vy_vs];
end


% FUNCTION: MVDR_LOCALISE_PRESSURE
% Standard MVDR scan over a 3D candidate grid using a pressure-only steering
% vector v = exp(-1i*k*R). Returns the position of the peak response.
function est = mvdr_localise_pressure(R, mic_xyz, grid_xyz, k, loading)
    [eigvecs, eigvals] = eig(R);
    eigvals = real(diag(eigvals));
    eigvals = max(eigvals, 1e-12);
    n_grid = size(grid_xyz, 2);
    out = zeros(n_grid, 1);
    for n = 1:n_grid
        l_pt = grid_xyz(:, n);
        r_lm = sqrt(sum((mic_xyz - l_pt).^2, 1)).';
        v = exp(-1i * k * r_lm);
        cbf = real(v' * R * v);
        lam_load = max(loading * cbf, loading * 1e-12);
        rxv = eigvecs * ((1.0 ./ (eigvals + lam_load)) .* (eigvecs' * v));
        denom = v' * rxv;
        if abs(denom) < 1e-12
            out(n) = 0;
        else
            w = rxv / denom;
            out(n) = abs(w' * R * w);
        end
    end
    [~, idx] = max(out);
    est = grid_xyz(:, idx);
end


% FUNCTION: MVDR_LOCALISE_AVS
% MVDR scan for AVS array with steering [P, u_x*P, u_y*P] per VS cluster
% (matches rho*c-scaled velocity convention used in vs_outputs).
% u = (mic - src) / |mic - src|, so V_x and V_y are the projection of
% the unit direction-of-propagation onto x and y axes (far-field convention).
function est = mvdr_localise_avs(R_avs, vs_centres, grid_xyz, k, loading)
    [eigvecs, eigvals] = eig(R_avs);
    eigvals = real(diag(eigvals));
    eigvals = max(eigvals, 1e-12);
    N_v = size(vs_centres, 2);
    n_grid = size(grid_xyz, 2);
    out = zeros(n_grid, 1);
    for n = 1:n_grid
        l_pt = grid_xyz(:, n);
        v_full = zeros(3 * N_v, 1);
        for vv = 1:N_v
            r_vec = vs_centres(:, vv) - l_pt;
            R = norm(r_vec);
            u = r_vec / R;
            P = exp(-1i * k * R);
            v_full(3*(vv-1) + (1:3)) = [P; u(1) * P; u(2) * P];
        end
        cbf = real(v_full' * R_avs * v_full);
        lam_load = max(loading * cbf, loading * 1e-12);
        rxv = eigvecs * ((1.0 ./ (eigvals + lam_load)) .* (eigvecs' * v_full));
        denom = v_full' * rxv;
        if abs(denom) < 1e-12
            out(n) = 0;
        else
            w = rxv / denom;
            out(n) = abs(w' * R_avs * w);
        end
    end
    [~, idx] = max(out);
    est = grid_xyz(:, idx);
end


% FUNCTION: FILL_BAND
% Shaded region between p25 and p75 percentile lines.
function fill_band(x, lo, hi, color, alpha)
    x = x(:).';
    xx = [x, fliplr(x)];
    yy = [lo(:).', fliplr(hi(:).')];
    fill(xx, yy, color, 'FaceAlpha', alpha, 'EdgeColor', 'none', ...
        'HandleVisibility', 'off');
end
