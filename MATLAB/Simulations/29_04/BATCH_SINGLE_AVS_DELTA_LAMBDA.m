clc; clear all; close all;

% Batch frequency × bearing sweep — Single Acoustic Vector Sensor MVDR
% Confirms the operational δ/λ band for the V₃ Spherical steering vector
% formulation, and characterises the FD velocity approximation error
% across the same (δ/λ, bearing) space.
%
% Sweep: frequencies × bearings at fixed δ and fixed source distance.
% Bearings are measured anticlockwise from the −x axis (so 0° = source
% along −x, matching the existing single-AVS test geometry).
%
% For each (f, bearing) point this script:
%   (1) Generates a synthetic signal via SUPER_VA_OUTPUT_BATCH.m
%   (2) Computes the FD velocity at the array centre and compares to the
%       analytical monopole particle velocity at (0,0,0) — diagnostic for
%       FD linearity breakdown (low r/λ wavefront curvature; high δ/λ
%       phase non-linearity)
%   (3) Runs MVDR beamforming, extracts radial error, angular error, and
%       beam pattern metrics
%   (4) Saves all results to CSV and .mat
%
% Output: line plots vs δ/λ (one line per bearing) for radial error,
% angular error, FD velocity error, beamwidth, and DI. Frequency points
% connected by lines; markers at each test frequency.
%
% NOTE on noise level:
% This script enforces noise_level = 0 in the generated signal — the
% goal is to map the *formulation* limits (FD breakdown, aliasing,
% near-field model) without contamination from SNR effects.
%
% Written by L Marshall 04/05/2026 — adapted from
% BATCH_FREQ_SWEEP_AVA_MVDR for single-AVS δ/λ characterisation.
%
% Patches (04/05/2026 evening):
%   - Angular error returns NaN when the cost function is saturated
%     (radial_error > 0.5 × margin_fixed). Previous version returned
%     numerically-zero angular error on principal-axis bearings due to
%     parabolic refinement landing on a grid cell with the exact correct
%     bearing — purely a coincidence of grid alignment, not a real metric
%   - FD velocity error normalisation fixed. Analytical reference is now
%     computed in the same normalised signal domain as the FD estimate,
%     by applying the same 1/sqrt(P_ref) scaling that the OUTPUT script
%     applies to the generated signal. Previous version compared
%     normalised FD against unnormalised analytical, giving meaningless
%     magnitude ratios of order 10–500
%   - Grating lobe detection added. For each beam pattern, counts peaks
%     within a configurable threshold of the mainlobe (default 3 dB) and
%     separated from it by more than the −3 dB beamwidth. Records count,
%     boolean presence, peak angle/separation/level. The first δ/λ at
%     which num_lobes ≥ 2 is the operational upper-band edge for
%     direction-finding without aliased ghosts
%   - 60° bearing dropped — geometrically equivalent to 30° by 4-fold
%     symmetry of the square mic layout under reflection across diagonals.
%     Previous run confirmed metrics match to 4 decimal places. Re-add by
%     editing batch_test_bearings_deg if the real array breaks symmetry
%
% Patches (continued):
%   - 3 dB range bracket along the source-bearing ray added as a per-test
%     diagnostic. Implementation matches Test A's `analyse_ray_peak`
%     conventions for direct comparability — same saturation criterion
%     (response within 3 dB of peak at BOTH ray endpoints), same parabolic
%     peak refinement, same 200-point ray sampling. Ray spans the same
%     physical range as the 2D search grid ([2δ, source_distance + margin])
%     so the bracket measurement is consistent with what the 2D MVDR scan
%     actually covers — bracket-saturation and 2D-scan-saturation agree
%     on what "saturated" means. Saturated points have radial AND angular
%     errors set to NaN. Ray response figures use Test A's plot style
%     (col_spherical brown-orange, LaTeX axis labels, southeast legend)
%     and saved at first/mid/last frequency per bearing, with the x-axis
%     in r/λ units to remain visually comparable to Test A figures.

% Test C and D


%% BATCH TEST PARAMETERS %%

% Test frequencies (Hz) — span below and above expected band edges to
% capture the failure regions on both sides
batch_test_frequencies = [630 800 1000 1250 1600 2000 2500 3150 4000];
batch_num_freqs = length(batch_test_frequencies);

% Test bearings (deg, measured anticlockwise from −x axis)
%   0°  → source along −x (broadside, sees δ side spacing)
%   45° → source along diagonal (sees δ√2 spacing — worst-case aliasing)
% NOTE: 60° bearing dropped from default sweep — geometrically equivalent
% to 30° by reflection across the y = −x diagonal (4-fold symmetry of the
% square BL/BR/TR/TL layout). Re-add to the array below if you need to
% probe asymmetry from real-array manufacturing tolerances.
batch_test_bearings_deg = [0, 15, 30, 45];
batch_num_bearings = length(batch_test_bearings_deg);

batch_num_tests = batch_num_freqs * batch_num_bearings;

% Fixed source distance (m)
source_distance_m = 0.4;

% Speed of sound, density
c_0 = 340;
rho_0 = 1.02;

% Array parameters
delta_fixed = 0.06; %MEMS colocation spacing (m)
N_a = 1; %single AVS — overrides anything in OUTPUT script
N_v = 1;
d_y = 0.1;

% Diagonal aperture — relevant Nyquist length for off-axis bearings
delta_diagonal = delta_fixed * sqrt(2);
f_nyq_side = c_0 / (2 * delta_fixed);
f_nyq_diag = c_0 / (2 * delta_diagonal);

% Processing parameters
overlap = 0.5;
loading = 1e-4;
analysis_halfwidth_hz = 10;
target_binwidth_hz = 1;

% Grating lobe detection — a peak counts as a grating lobe if:
%   (a) it is a local maximum in the polar beam pattern
%   (b) its level is within `grating_lobe_threshold_db` of the mainlobe
%   (c) it is separated from the mainlobe by more than the mainlobe's
%       −3 dB half-beamwidth (so we are not just counting points within
%       the mainlobe lobe-shoulder)
% The first δ/λ at which num_lobes ≥ 2 is the operational upper-band
% edge for direction-finding without aliased ghosts. Standard convention
% in beamforming literature is 3 dB ("another peak as strong as the
% mainlobe"); 6 dB is a more lenient criterion.
grating_lobe_threshold_db = 3.0;

% 3 dB range bracket parameters
% 200-point ray sampling matches Test A. The ray spans the same physical
% range as the 2D search grid: [2δ, source_distance + margin].
ray_num_points = 200;

% Grid parameters
% For confirming sub-grid radial errors at high δ/λ, prefer a fine fixed
% grid. Rule of thumb: grid_resolution ≤ 0.005 m gives ~22 pts/λ at
% 3 kHz, enough for parabolic refinement to be meaningful.
grid_pts_per_lambda = 20;
grid_resolution_fixed = 0.005; %m, [] for λ-scaled
margin_fixed = 0.5;

% Force noise-free signal generation — formulation envelope test
batch_test_noise_level = 0.0;

% Create results folder
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
results_folder = sprintf('single_avs_dl_sweep_%s', timestamp);
mkdir(results_folder);

fprintf('\n<strong>SINGLE-AVS δ/λ × BEARING SWEEP</strong>\n');
fprintf('Frequencies (%d): ', batch_num_freqs);
fprintf('%.0f  ', batch_test_frequencies);
fprintf('Hz\n');
fprintf('Bearings (%d): ', batch_num_bearings);
fprintf('%.0f°  ', batch_test_bearings_deg);
fprintf('(from −x axis)\n');
fprintf('Source distance: %.3f m\n', source_distance_m);
fprintf('δ (side): %.4f m | δ_diag: %.4f m\n', delta_fixed, delta_diagonal);
fprintf('f_Nyquist side: %.0f Hz | f_Nyquist diagonal: %.0f Hz\n', f_nyq_side, f_nyq_diag);
fprintf('δ/λ range: [%.3f, %.3f]\n', ...
    delta_fixed * min(batch_test_frequencies)/c_0, ...
    delta_fixed * max(batch_test_frequencies)/c_0);
fprintf('Noise level: %.3f (formulation envelope test)\n', batch_test_noise_level);
if isempty(grid_resolution_fixed)
    fprintf('Grid: %d pts/λ (frequency-scaled) | Margin: %.2f m\n', ...
        grid_pts_per_lambda, margin_fixed);
else
    fprintf('Grid: %.4f m fixed | Margin: %.2f m\n', grid_resolution_fixed, margin_fixed);
end
fprintf('Total tests: %d\n', batch_num_tests);
fprintf('Results folder: %s\n\n', results_folder);


%% INITIALISE RESULTS STORAGE %%
% Stored as flat vectors of length batch_num_tests; bearing and δ/λ are
% kept as separate columns so the post-processing can group either way.

batch_res_test_idx = zeros(batch_num_tests, 1);
batch_res_freq_hz = zeros(batch_num_tests, 1);
batch_res_lambda_m = zeros(batch_num_tests, 1);
batch_res_dl_ratio = zeros(batch_num_tests, 1);
batch_res_dl_diag_ratio = zeros(batch_num_tests, 1);
batch_res_bearing_deg = zeros(batch_num_tests, 1);
batch_res_source_x = zeros(batch_num_tests, 1);
batch_res_source_y = zeros(batch_num_tests, 1);
batch_res_grid_resolution = zeros(batch_num_tests, 1);

% Localisation metrics
batch_res_radial_error = zeros(batch_num_tests, 1);
batch_res_angular_error = zeros(batch_num_tests, 1);
batch_res_peak_db = zeros(batch_num_tests, 1);

% 3 dB range bracket along the source-bearing ray. Ray spans the same
% physical range as the 2D search grid ([2δ, source_distance + margin]),
% so the bracket measurement matches what the 2D MVDR scan covers.
% Two saturation criteria stored separately:
%   *_testA: Test A style — response within 3 dB of peak at BOTH
%            endpoints of the ray. Catches "flat cost function" failure
%   *_peakOff: peak displacement — ray peak landed > 0.5 × margin from
%              the true source. Catches "peak migrated to boundary"
% The canonical `Saturated` column is OR of the two.
batch_res_bracket_m = zeros(batch_num_tests, 1);
batch_res_bracket_rl = zeros(batch_num_tests, 1);
batch_res_saturated = false(batch_num_tests, 1);
batch_res_saturated_testA = false(batch_num_tests, 1);
batch_res_saturated_peakOff = false(batch_num_tests, 1);
batch_res_ray_peak_r = zeros(batch_num_tests, 1);   %peak position along ray (m)
batch_res_peak_displacement = zeros(batch_num_tests, 1);  %|ray_peak − true r| (m)

% Beam pattern metrics
batch_res_beamwidth_deg = zeros(batch_num_tests, 1);
batch_res_sidelobe_db = zeros(batch_num_tests, 1);
batch_res_DI_db = zeros(batch_num_tests, 1);

% Grating lobe metrics
batch_res_num_lobes = zeros(batch_num_tests, 1);
batch_res_grating_present = false(batch_num_tests, 1);
batch_res_grating_angle_deg = nan(batch_num_tests, 1);
batch_res_grating_separation_deg = nan(batch_num_tests, 1);
batch_res_grating_level_db = nan(batch_num_tests, 1);

% FD velocity error metrics — analytical-field diagnostic (Test B style)
% Two reported quantities:
%   _rel — |v_fd − v_true| / |v_true| (relative norm, scale-invariant)
%   _ang_deg — angle between FD-estimated and true velocity vectors (deg)
batch_res_fd_vel_err_rel = zeros(batch_num_tests, 1);
batch_res_fd_vel_err_ang_deg = zeros(batch_num_tests, 1);


%% RUN BATCH TESTS %%

test_counter = 0;

for bearing_idx = 1:batch_num_bearings
    bearing_deg = batch_test_bearings_deg(bearing_idx);
    bearing_rad = deg2rad(bearing_deg);

    % Bearings measured anticlockwise from −x axis. Convention:
    %   x = −R cos(bearing), y = −R sin(bearing)
    % This puts 0° at (−R, 0), 45° at (−R/√2, −R/√2), 90° at (0, −R).
    source_x = -source_distance_m * cos(bearing_rad);
    source_y = -source_distance_m * sin(bearing_rad);

    fprintf('\n##############################################\n');
    fprintf('<strong>BEARING %d/%d: %.0f° → source (%.3f, %.3f) m</strong>\n', ...
        bearing_idx, batch_num_bearings, bearing_deg, source_x, source_y);
    fprintf('##############################################\n');

    for freq_idx = 1:batch_num_freqs

        test_counter = test_counter + 1;
        batch_test_freq = batch_test_frequencies(freq_idx);
        lambda = c_0 / batch_test_freq;
        dl_ratio = delta_fixed / lambda;
        dl_diag_ratio = delta_diagonal / lambda;
        rl_ratio = source_distance_m / lambda;

        fprintf('\n========================================\n');
        fprintf('<strong>TEST %d/%d (bearing %.0f°, %.0f Hz)</strong>\n', ...
            test_counter, batch_num_tests, bearing_deg, batch_test_freq);
        fprintf('  λ = %.4f m | δ/λ = %.3f | δ_diag/λ = %.3f | r/λ = %.2f\n', ...
            lambda, dl_ratio, dl_diag_ratio, rl_ratio);
        fprintf('========================================\n');

        % Grid step
        if isempty(grid_resolution_fixed)
            grid_resolution = lambda / grid_pts_per_lambda;
        else
            grid_resolution = grid_resolution_fixed;
        end

        f_band_lo = batch_test_freq - analysis_halfwidth_hz;
        f_band_hi = batch_test_freq + analysis_halfwidth_hz;

        % Store sweep coordinates
        batch_res_test_idx(test_counter) = test_counter;
        batch_res_freq_hz(test_counter) = batch_test_freq;
        batch_res_lambda_m(test_counter) = lambda;
        batch_res_dl_ratio(test_counter) = dl_ratio;
        batch_res_dl_diag_ratio(test_counter) = dl_diag_ratio;
        batch_res_bearing_deg(test_counter) = bearing_deg;
        batch_res_source_x(test_counter) = source_x;
        batch_res_source_y(test_counter) = source_y;
        batch_res_grid_resolution(test_counter) = grid_resolution;


        % STEP 1: CLEAR PREVIOUS RUN VARIABLES %

        vars_to_clear = {'source_positions', 'source_frequencies', ...
            'source_amplitudes', 'all_results'};
        for v = 1:length(vars_to_clear)
            if exist(vars_to_clear{v}, 'var')
                clear(vars_to_clear{v});
            end
        end


        % STEP 2: RUN SIGNAL GENERATION SCRIPT %
        % NB: SUPER_VA_OUTPUT_BATCH.m must have N_a = 1 set internally
        % for single-AVS output (4 mic columns). If your version of that
        % script defaults to N_a = 2, copy it to a single-AVS variant
        % first or edit the N_a default.

        fprintf('\nGenerating signal (noise = %.3f)...\n', batch_test_noise_level);

        batch_test_delta = delta_fixed;
        batch_test_source_x = source_x;
        batch_test_source_y = source_y;
        batch_test_distance = source_distance_m;
        batch_test_distance_lambda = rl_ratio;
        batch_csv_name = fullfile(results_folder, ...
            sprintf('signal_b%02d_f%02d_%.0fhz.csv', ...
            bearing_idx, freq_idx, batch_test_freq));

        % noise_level override is read by SUPER_VA_OUTPUT_BATCH if the
        % batch_test_noise_level variable exists when it runs

        run('SUPER_VA_OUTPUT_BATCH.m');
        fprintf('  Signal generation complete\n');


        % STEP 3: LOAD SIGNAL AND CHECK FORMAT %

        data_table = readtable(csv_filename);
        data = table2array(data_table);
        time = data(:, 1);
        N = size(data, 1);
        N_ma = N_v * 4;

        N_m_total = size(data, 2) - 1;
        if N_m_total ~= N_a * N_ma
            warning('CSV has %d mic cols, expected %d (N_a=%d, N_v=%d). OUTPUT script may have wrong N_a.', ...
                N_m_total, N_a*N_ma, N_a, N_v);
        end

        % Single AVS — the 4 mic columns
        tx_vs = data(:, 2:1+N_ma).';

        fprintf('  Loaded %d samples × %d mics\n', N, N_ma);


        % STEP 4: ARRAY GEOMETRY %

        A1_coord = [0; 0; 0];
        mic_offsets = delta_fixed / 2 * [-1, 1, 1, -1; -1, -1, 1, 1; 0, 0, 0, 0];

        vs_centres = zeros(3, N_v);
        mic_positions = zeros(3, N_a * N_ma);
        for vs = 1:N_v
            vs_centres(:, vs) = A1_coord + [0; (vs - 1) * d_y; 0];
            idx_mics = (vs - 1) * 4 + (1:4);
            mic_positions(:, idx_mics) = vs_centres(:, vs) + mic_offsets;
        end


        % STEP 5: FFT SETUP %

        d_t = time(2) - time(1);
        F_s = 1 / d_t;

        cycles_per_window = round(batch_test_freq / target_binwidth_hz);
        size_fft = round(cycles_per_window * F_s / batch_test_freq);
        actual_binwidth = F_s / size_fft;

        if size_fft > N / 4
            size_fft = 2^floor(log2(N / 4));
            actual_binwidth = F_s / size_fft;
        end

        fft_vec = F_s * (0:(size_fft - 1)) / size_fft;
        bin_mask = (fft_vec >= f_band_lo) & (fft_vec <= f_band_hi);
        bin_index = find(bin_mask);
        bin_freqs = fft_vec(bin_index);
        num_bins = length(bin_index);

        if num_bins == 0
            [~, nearest] = min(abs(fft_vec(1:floor(size_fft/2)) - batch_test_freq));
            bin_index = nearest;
            bin_freqs = fft_vec(nearest);
            num_bins = 1;
        end

        fprintf('  size_fft = %d | bin width = %.3f Hz | %d bins\n', ...
            size_fft, actual_binwidth, num_bins);


        % STEP 6: SNAPSHOTS AND CSM %

        window = hanning(size_fft)';
        snapshots_vs = make_snapshots(tx_vs, size_fft, overlap, window);
        r_vs = create_vs_csm(snapshots_vs, bin_index, delta_fixed, ...
            bin_freqs, rho_0, N_v, num_bins, c_0);

        fprintf('  CSM: %d × %d × %d\n', size(r_vs, 1), size(r_vs, 2), num_bins);


        % STEP 7: FD VELOCITY ERROR DIAGNOSTIC %
        % Analytical-field diagnostic, matching Test B's
        % compute_velocity_fd_error.  Synthesises the analytical monopole
        % pressure phasor at each mic position, runs the same
        % least-squares FD fit that vs_outputs uses internally, and
        % compares against the analytical Euler-equation velocity at the
        % AVS centre. Both quantities live in the same physical units —
        % no signal normalisation chain to debug.
        %
        % This isolates the geometric question (does the linear-gradient
        % assumption hold across the array?) from the simulated-signal
        % chain. Numbers are directly comparable to Test B values.
        [fd_rel_err, fd_ang_err_deg] = compute_velocity_fd_error( ...
            [source_x; source_y; 0], mic_positions, vs_centres, ...
            batch_test_freq, c_0, rho_0, N_v);

        batch_res_fd_vel_err_rel(test_counter) = fd_rel_err;
        batch_res_fd_vel_err_ang_deg(test_counter) = fd_ang_err_deg;

        fprintf('\n<strong>FD velocity error (analytical, f = %.0f Hz):</strong>\n', batch_test_freq);
        fprintf('  |v_fd − v_true|/|v_true| = %.4f\n', fd_rel_err);
        fprintf('  Angular error (FD vs true vector): %.3f°\n', fd_ang_err_deg);


        % STEP 8: SEARCH GRID %

        system_centre_y = 0; %single-AVS — array centre is origin
        [x_scan, y_scan, X_grid, Y_grid, candidate_points, grid_res] = ...
            build_search_grid(batch_test_freq, c_0, source_distance_m, ...
            system_centre_y, grid_pts_per_lambda, margin_fixed, grid_resolution_fixed);

        fprintf('  Search grid: %d × %d points, resolution %.4f m\n', ...
            length(x_scan), length(y_scan), grid_res);


        % STEP 9: MVDR BEAMFORMING %

        fprintf('\n<strong>Running MVDR...</strong>\n');
        source_positions = [source_x, source_y, 0];

        response_db = mvdr_vs_beamforming(r_vs, vs_centres, candidate_points, ...
            bin_freqs, c_0, rho_0, loading, num_bins, N_v);
        response_db = response_db - max(response_db);

        fprintf('  Response range: [%.2f, %.2f] dB\n', ...
            min(response_db), max(response_db));


        % STEP 10: LOCALISATION METRICS + 3 dB BRACKET %

        grid_response = reshape(response_db, length(y_scan), length(x_scan));
        [est_x, est_y, peak_db] = refine_peak_2d(grid_response, x_scan, y_scan);
        est_pos = [est_x, est_y];

        true_pos = [source_x, source_y];
        radial_err = norm(true_pos - est_pos);

        % Compute the 3 dB range bracket along the source-bearing ray.
        % Ray spans the same physical range as the 2D search grid so the
        % bracket measurement is consistent with what the 2D MVDR scan
        % actually covers — bracket-saturation and 2D-saturation agree.
        max_extent = source_distance_m + margin_fixed;
        bearing_2d = [cos(bearing_rad); sin(bearing_rad)];
        % Ray direction from origin toward the source
        ray_dir = -bearing_2d; % source is at (-R cos(b), -R sin(b))
        ray_lo = 2 * delta_fixed;
        ray_hi = max_extent;

        [bracket_m, bracket_rl, ray_peak_r, ray_response_db, ray_r_vals, saturated_testA] = ...
            compute_range_bracket(r_vs, vs_centres, ray_dir, ...
            bin_freqs, c_0, rho_0, loading, num_bins, N_v, ...
            ray_lo, ray_hi, lambda, ray_num_points);

        % Saturation detection — combine two failure modes:
        %   (1) Test A criterion: response within 3 dB of peak at BOTH ray
        %       endpoints. Catches the "cost function flat across whole
        %       search range" mode (relevant when r is varied and the
        %       cost function never localises).
        %   (2) Peak-displacement criterion: the ray peak landed near the
        %       search boundary rather than near the true source. Catches
        %       the "peak migrated to boundary" mode that fails Test A's
        %       criterion (when the cost function has a peak but at the
        %       wrong place — common when r is fixed and δ/λ is varied).
        % This sweep needs (2); Test A originally needed (1). Keeping both
        % as separate columns lets either test be scored on its own terms.
        peak_displacement = abs(ray_peak_r - source_distance_m);
        saturated_peakOff = peak_displacement > 0.5 * margin_fixed;

        saturation_flag = saturated_testA || saturated_peakOff;

        % Angular error: only meaningful when the cost function has a
        % well-localised peak in 2D. NaN when saturated.
        ref_pos = [0, 0]; %array centre
        if saturation_flag
            angular_err_deg = NaN;
        else
            true_angle = atan2(true_pos(2) - ref_pos(2), true_pos(1) - ref_pos(1));
            est_angle  = atan2(est_pos(2) - ref_pos(2), est_pos(1) - ref_pos(1));
            angular_err_deg = rad2deg(abs(true_angle - est_angle));
            if angular_err_deg > 180
                angular_err_deg = 360 - angular_err_deg;
            end
        end

        % Also NaN out the radial error itself when saturated — the
        % radial error of a saturated peak is meaningless (peak position
        % is determined by numerical noise on a flat plateau)
        if saturation_flag
            radial_err_reported = NaN;
        else
            radial_err_reported = radial_err;
        end

        batch_res_radial_error(test_counter) = radial_err_reported;
        batch_res_angular_error(test_counter) = angular_err_deg;
        batch_res_peak_db(test_counter) = peak_db;
        batch_res_bracket_m(test_counter) = bracket_m;
        batch_res_bracket_rl(test_counter) = bracket_rl;
        batch_res_saturated(test_counter) = saturation_flag;
        batch_res_saturated_testA(test_counter) = saturated_testA;
        batch_res_saturated_peakOff(test_counter) = saturated_peakOff;
        batch_res_ray_peak_r(test_counter) = ray_peak_r;
        batch_res_peak_displacement(test_counter) = peak_displacement;

        fprintf('  True: (%.4f, %.4f) | Est: (%.4f, %.4f)\n', ...
            true_pos, est_pos);
        fprintf('  3 dB bracket: %.4f m (%.3f r/λ)\n', bracket_m, bracket_rl);
        fprintf('  Ray peak: r = %.4f m (true = %.4f m, displacement = %.4f m)\n', ...
            ray_peak_r, source_distance_m, peak_displacement);
        if saturation_flag
            sat_reasons = {};
            if saturated_testA, sat_reasons{end+1} = 'flat-CF'; end
            if saturated_peakOff, sat_reasons{end+1} = 'peak-off-bdry'; end
            fprintf('  SATURATED [%s] — radial and angular errors set NaN\n', ...
                strjoin(sat_reasons, ', '));
        else
            fprintf('  Radial err: %.4f m | Angular err: %.3f°\n', ...
                radial_err, angular_err_deg);
        end

        % Save ray response figure on selected test points (same trigger
        % as the 2D scan figures — first, mid, last frequency per bearing).
        % Style matches Test A's plot_ray_response: 520x360 figure, brown-
        % orange palette (col_spherical for the single V_3 line in this
        % single-formulation test), LaTeX-typeset axis labels, southeast
        % legend, dashed black true-distance line. X-axis in r/λ keeps
        % the visual comparison to Test A even though the ray range is
        % set in metres (matched to the 2D search grid).
        save_ray_response = (freq_idx == 1) || ...
                            (freq_idx == round(batch_num_freqs/2)) || ...
                            (freq_idx == batch_num_freqs);
        if save_ray_response
            col_spherical = [0.550, 0.300, 0.080]; %dark orange-brown
            fn_ray = 'Times New Roman';

            fig_ray = figure('Color', 'w', 'Position', [100 100 520 360]);
            plot(ray_r_vals / lambda, ray_response_db, '-', ...
                'Color', col_spherical, 'LineWidth', 1.8, ...
                'DisplayName', 'Spherical $\mathbf{V}_3$');
            hold on;
            xline(source_distance_m / lambda, 'k--', 'LineWidth', 1.2, ...
                'DisplayName', 'True distance', 'HandleVisibility', 'on');
            yline(-3, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0, ...
                'HandleVisibility', 'off');
            hold off;
            xlabel('Distance along ray ($r/\lambda$)', ...
                'Interpreter', 'latex', 'FontSize', 14);
            ylabel('Normalised response (dB)', 'FontName', fn_ray, 'FontSize', 14);
            legend('Location', 'southeast', 'FontName', fn_ray, 'FontSize', 11, ...
                'Interpreter', 'latex', 'Box', 'on');
            set(gca, 'FontName', fn_ray, 'FontSize', 12, ...
                'LineWidth', 1.0, 'Box', 'on', 'YGrid', 'on', 'XGrid', 'off');
            title(sprintf('Ray response — bearing %.0f°, f = %.0f Hz, $\\delta/\\lambda$ = %.3f', ...
                bearing_deg, batch_test_freq, dl_ratio), ...
                'FontName', fn_ray, 'FontSize', 12, 'FontWeight', 'normal', ...
                'Interpreter', 'latex');
            exportgraphics(gcf, fullfile(results_folder, ...
                sprintf('ray_b%02d_f%02d_%.0fhz.png', ...
                bearing_idx, freq_idx, batch_test_freq)), 'ContentType', 'image');
            exportgraphics(gcf, fullfile(results_folder, ...
                sprintf('ray_b%02d_f%02d_%.0fhz.pdf', ...
                bearing_idx, freq_idx, batch_test_freq)), 'ContentType', 'vector');
            close(gcf);
        end


        % STEP 11: 2D SCAN PLOT (ON SELECTED TEST POINTS ONLY) %
        % Saving all 36 spatial maps would clutter the results folder.
        % Save only at the lowest, mid, and highest f at each bearing.

        save_2d_scan = (freq_idx == 1) || ...
                       (freq_idx == round(batch_num_freqs/2)) || ...
                       (freq_idx == batch_num_freqs);

        if save_2d_scan
            plot_2dscan_single(grid_response, x_scan, y_scan, ...
                source_positions, est_pos, vs_centres, bearing_deg, batch_test_freq);
            exportgraphics(gcf, fullfile(results_folder, ...
                sprintf('scan_b%02d_f%02d_%.0fhz.png', ...
                bearing_idx, freq_idx, batch_test_freq)), 'ContentType', 'image');
            close(gcf);
        end


        % STEP 12: BEAM PATTERN METRICS %

        fprintf('\n<strong>Beam pattern metrics...</strong>\n');
        [beam_metrics, fig_handle] = compute_nearfield_beam_pattern(r_vs, ...
            vs_centres, mean(vs_centres(1:2,:), 2), source_positions(:,1:2), ...
            bin_freqs, c_0, rho_0, loading, num_bins, N_v, ...
            source_distance_m, ...
            sprintf('Beam pattern (b=%.0f°, f=%.0f Hz)', bearing_deg, batch_test_freq), ...
            grating_lobe_threshold_db);

        if save_2d_scan
            exportgraphics(fig_handle, fullfile(results_folder, ...
                sprintf('beam_b%02d_f%02d_%.0fhz.png', ...
                bearing_idx, freq_idx, batch_test_freq)), 'ContentType', 'image');
        end
        close(fig_handle);

        batch_res_beamwidth_deg(test_counter) = beam_metrics.beamwidth_3dB;
        batch_res_sidelobe_db(test_counter) = beam_metrics.sidelobe_level;
        batch_res_DI_db(test_counter) = beam_metrics.directivity_index;

        % Grating lobe metrics
        batch_res_num_lobes(test_counter) = beam_metrics.num_lobes;
        batch_res_grating_present(test_counter) = beam_metrics.grating_lobe_present;
        batch_res_grating_angle_deg(test_counter) = beam_metrics.grating_lobe_angle_deg;
        batch_res_grating_separation_deg(test_counter) = beam_metrics.grating_lobe_separation_deg;
        batch_res_grating_level_db(test_counter) = beam_metrics.grating_lobe_level_db;

    end % freq loop
end % bearing loop


%% SAVE RESULTS %%

fprintf('\n========================================\n');
fprintf('<strong>SAVING RESULTS</strong>\n');
fprintf('========================================\n');

results_table = table( ...
    batch_res_test_idx, batch_res_bearing_deg, batch_res_freq_hz, ...
    batch_res_lambda_m, batch_res_dl_ratio, batch_res_dl_diag_ratio, ...
    batch_res_source_x, batch_res_source_y, batch_res_grid_resolution, ...
    batch_res_radial_error, batch_res_angular_error, batch_res_peak_db, ...
    batch_res_bracket_m, batch_res_bracket_rl, ...
    batch_res_saturated, batch_res_saturated_testA, batch_res_saturated_peakOff, ...
    batch_res_ray_peak_r, batch_res_peak_displacement, ...
    batch_res_beamwidth_deg, batch_res_sidelobe_db, batch_res_DI_db, ...
    batch_res_num_lobes, batch_res_grating_present, ...
    batch_res_grating_angle_deg, batch_res_grating_separation_deg, ...
    batch_res_grating_level_db, ...
    batch_res_fd_vel_err_rel, batch_res_fd_vel_err_ang_deg, ...
    'VariableNames', {'TestIdx', 'Bearing_deg', 'Freq_Hz', ...
    'Lambda_m', 'DeltaLambdaSide', 'DeltaLambdaDiag', ...
    'Source_x_m', 'Source_y_m', 'Grid_res_m', ...
    'Radial_err_m', 'Angular_err_deg', 'Peak_dB', ...
    'Bracket_3dB_m', 'Bracket_3dB_rlambda', ...
    'Saturated', 'Saturated_testA_flatCF', 'Saturated_peakOffBoundary', ...
    'Ray_peak_r_m', 'Peak_displacement_m', ...
    'Beamwidth_3dB_deg', 'Sidelobe_dB', 'DI_dB', ...
    'Num_lobes', 'Grating_present', ...
    'Grating_angle_deg', 'Grating_separation_deg', 'Grating_level_dB', ...
    'FDvel_rel_err', 'FDvel_ang_err_deg'});

disp(results_table);

writetable(results_table, fullfile(results_folder, 'sweep_results.csv'));
save(fullfile(results_folder, 'sweep_results.mat'), 'results_table', ...
    'batch_test_frequencies', 'batch_test_bearings_deg', ...
    'source_distance_m', 'delta_fixed', 'delta_diagonal', ...
    'f_nyq_side', 'f_nyq_diag', 'c_0', 'rho_0');

fprintf('\nSaved: %s\n', fullfile(results_folder, 'sweep_results.csv'));


%% PUBLICATION FIGURES — VS δ/λ %%

fn = 'Times New Roman';
fs_ax = 12; fs_lab = 14; fs_leg = 11;
lw = 1.6; ms = 7;

bearing_colours = lines(batch_num_bearings);
bearing_markers = {'o', 's', '^', 'd', 'v', '>'};

% Sort each bearing's data by δ/λ for clean line plots
bearing_data = cell(batch_num_bearings, 1);
for bi = 1:batch_num_bearings
    mask = batch_res_bearing_deg == batch_test_bearings_deg(bi);
    [dl_sorted, sort_idx] = sort(batch_res_dl_ratio(mask));
    bd = struct();
    bd.dl = dl_sorted;
    bd.dl_diag = batch_res_dl_diag_ratio(mask); bd.dl_diag = bd.dl_diag(sort_idx);
    bd.freq = batch_res_freq_hz(mask); bd.freq = bd.freq(sort_idx);
    bd.radial = batch_res_radial_error(mask); bd.radial = bd.radial(sort_idx);
    bd.angular = batch_res_angular_error(mask); bd.angular = bd.angular(sort_idx);
    bd.beamwidth = batch_res_beamwidth_deg(mask); bd.beamwidth = bd.beamwidth(sort_idx);
    bd.di = batch_res_DI_db(mask); bd.di = bd.di(sort_idx);
    bd.sidelobe = batch_res_sidelobe_db(mask); bd.sidelobe = bd.sidelobe(sort_idx);
    bd.fd_err = batch_res_fd_vel_err_rel(mask); bd.fd_err = bd.fd_err(sort_idx);
    bd.fd_ang = batch_res_fd_vel_err_ang_deg(mask); bd.fd_ang = bd.fd_ang(sort_idx);
    bd.num_lobes = batch_res_num_lobes(mask); bd.num_lobes = bd.num_lobes(sort_idx);
    bd.grating_present = batch_res_grating_present(mask); bd.grating_present = bd.grating_present(sort_idx);
    bd.grating_level = batch_res_grating_level_db(mask); bd.grating_level = bd.grating_level(sort_idx);
    bd.bracket_m = batch_res_bracket_m(mask); bd.bracket_m = bd.bracket_m(sort_idx);
    bd.bracket_rl = batch_res_bracket_rl(mask); bd.bracket_rl = bd.bracket_rl(sort_idx);
    bd.saturated = batch_res_saturated(mask); bd.saturated = bd.saturated(sort_idx);
    bearing_data{bi} = bd;
end

% Helper: plot a metric vs δ/λ with one line per bearing
plot_metric_vs_dl = @(metric_field, ylabel_str, use_log, save_name) ...
    do_plot_metric(bearing_data, batch_test_bearings_deg, bearing_colours, ...
    bearing_markers, metric_field, ylabel_str, use_log, ...
    fullfile(results_folder, save_name), fn, fs_ax, fs_lab, fs_leg, lw, ms, ...
    delta_fixed, c_0, f_nyq_side, f_nyq_diag);


% RADIAL ERROR
plot_metric_vs_dl('radial', 'Radial error (m)', true, 'radial_error_vs_dl.png');

% ANGULAR ERROR
plot_metric_vs_dl('angular', 'Angular error (deg)', true, 'angular_error_vs_dl.png');

% FD VELOCITY ERROR (relative)
plot_metric_vs_dl('fd_err', '|v_{FD} - v_{true}| / |v_{true}|', true, 'fd_velocity_err_vs_dl.png');

% FD VELOCITY ANGULAR ERROR
plot_metric_vs_dl('fd_ang', 'FD velocity vector angular error (deg)', false, 'fd_velocity_ang_err_vs_dl.png');

% BEAMWIDTH
plot_metric_vs_dl('beamwidth', '−3 dB beamwidth (deg)', false, 'beamwidth_vs_dl.png');

% DIRECTIVITY INDEX
plot_metric_vs_dl('di', 'Directivity index (dB)', false, 'DI_vs_dl.png');


% GRATING LOBE: 2-panel figure replicating the user's reference layout
% (top: number of lobes detected; bottom: grating lobe present yes/no)
figure('Color', 'w', 'Position', [100 100 700 500]);

subplot(2,1,1);
for bi = 1:batch_num_bearings
    bd = bearing_data{bi};
    plot(bd.dl, bd.num_lobes, '-', 'Color', bearing_colours(bi,:), ...
        'LineWidth', lw, 'Marker', bearing_markers{bi}, 'MarkerSize', ms, ...
        'MarkerFaceColor', bearing_colours(bi,:), ...
        'DisplayName', sprintf('%.0f°', batch_test_bearings_deg(bi)));
    hold on;
end
xline(delta_fixed * f_nyq_diag / c_0, '--k', 'HandleVisibility', 'off');
ylabel('Number of lobes', 'FontName', fn, 'FontSize', fs_lab);
title('Number of lobes detected', 'FontName', fn, 'FontWeight', 'normal');
legend('Location', 'best', 'FontName', fn, 'FontSize', fs_leg, 'Box', 'on');
set(gca, 'FontName', fn, 'FontSize', fs_ax, 'YGrid', 'on', 'Box', 'on', ...
    'YTick', 0:max(cellfun(@(b) max([b.num_lobes; 0]), bearing_data)));

subplot(2,1,2);
for bi = 1:batch_num_bearings
    bd = bearing_data{bi};
    stem(bd.dl, double(bd.grating_present), '-', 'Color', bearing_colours(bi,:), ...
        'LineWidth', lw, 'Marker', bearing_markers{bi}, 'MarkerSize', ms, ...
        'MarkerFaceColor', bearing_colours(bi,:), ...
        'DisplayName', sprintf('%.0f°', batch_test_bearings_deg(bi)));
    hold on;
end
xline(delta_fixed * f_nyq_diag / c_0, '--k', 'HandleVisibility', 'off');
xlabel('\delta / \lambda', 'FontName', fn, 'FontSize', fs_lab);
ylabel('Grating lobe present', 'FontName', fn, 'FontSize', fs_lab);
title('Grating lobe detection', 'FontName', fn, 'FontWeight', 'normal');
ylim([-0.1, 1.2]);
set(gca, 'FontName', fn, 'FontSize', fs_ax, 'Box', 'on', ...
    'YTick', [0, 1], 'YTickLabel', {'No', 'Yes'});

exportgraphics(gcf, fullfile(results_folder, 'grating_lobe_detection_vs_dl.png'), ...
    'ContentType', 'image');


% Strongest grating-lobe level (NaN where none present — gaps in line)
plot_metric_vs_dl('grating_level', 'Strongest grating lobe level (dB)', ...
    false, 'grating_level_vs_dl.png');


% 3 dB BRACKET vs δ/λ — primary saturation diagnostic. Ray spans the
% same physical range as the 2D search grid.
plot_metric_vs_dl('bracket_rl', '3 dB range bracket (r/\lambda)', ...
    true, 'bracket_vs_dl.png');


%% COMPOSITE 2x3 SUMMARY FIGURE %%

figure('Color', 'w', 'Position', [100 100 1300 700]);

% (a) Radial error
subplot(2,3,1);
for bi = 1:batch_num_bearings
    bd = bearing_data{bi};
    semilogy(bd.dl, bd.radial, '-', 'Color', bearing_colours(bi,:), ...
        'LineWidth', lw, 'Marker', bearing_markers{bi}, 'MarkerSize', ms-1, ...
        'MarkerFaceColor', bearing_colours(bi,:), ...
        'DisplayName', sprintf('%.0f°', batch_test_bearings_deg(bi)));
    hold on;
end
xline(delta_fixed * f_nyq_diag / c_0, '--k', 'HandleVisibility', 'off');
xlabel('\delta / \lambda', 'FontName', fn, 'FontSize', fs_lab);
ylabel('Radial error (m)', 'FontName', fn, 'FontSize', fs_lab);
legend('Location', 'best', 'FontName', fn, 'FontSize', fs_leg, 'Box', 'on');
set(gca, 'FontName', fn, 'FontSize', fs_ax, 'YGrid', 'on', 'Box', 'on');
title('(a) Radial error (NaN = saturated)', 'FontName', fn, 'FontWeight', 'normal');

% (b) Angular error
subplot(2,3,2);
for bi = 1:batch_num_bearings
    bd = bearing_data{bi};
    semilogy(bd.dl, max(bd.angular, 1e-3), '-', 'Color', bearing_colours(bi,:), ...
        'LineWidth', lw, 'Marker', bearing_markers{bi}, 'MarkerSize', ms-1, ...
        'MarkerFaceColor', bearing_colours(bi,:), ...
        'DisplayName', sprintf('%.0f°', batch_test_bearings_deg(bi)));
    hold on;
end
xline(delta_fixed * f_nyq_diag / c_0, '--k', 'HandleVisibility', 'off');
xlabel('\delta / \lambda', 'FontName', fn, 'FontSize', fs_lab);
ylabel('Angular error (deg)', 'FontName', fn, 'FontSize', fs_lab);
set(gca, 'FontName', fn, 'FontSize', fs_ax, 'YGrid', 'on', 'Box', 'on');
title('(b) Angular error', 'FontName', fn, 'FontWeight', 'normal');

% (c) 3 dB range bracket — primary saturation diagnostic
subplot(2,3,3);
for bi = 1:batch_num_bearings
    bd = bearing_data{bi};
    semilogy(bd.dl, bd.bracket_rl, '-', 'Color', bearing_colours(bi,:), ...
        'LineWidth', lw, 'Marker', bearing_markers{bi}, 'MarkerSize', ms-1, ...
        'MarkerFaceColor', bearing_colours(bi,:), ...
        'DisplayName', sprintf('%.0f°', batch_test_bearings_deg(bi)));
    hold on;
end
xline(delta_fixed * f_nyq_diag / c_0, '--k', 'HandleVisibility', 'off');
xlabel('\delta / \lambda', 'FontName', fn, 'FontSize', fs_lab);
ylabel('3 dB bracket (r/\lambda)', 'FontName', fn, 'FontSize', fs_lab);
set(gca, 'FontName', fn, 'FontSize', fs_ax, 'YGrid', 'on', 'Box', 'on');
title('(c) 3 dB range bracket', 'FontName', fn, 'FontWeight', 'normal');

% (d) FD velocity error
subplot(2,3,4);
for bi = 1:batch_num_bearings
    bd = bearing_data{bi};
    semilogy(bd.dl, bd.fd_err, '-', 'Color', bearing_colours(bi,:), ...
        'LineWidth', lw, 'Marker', bearing_markers{bi}, 'MarkerSize', ms-1, ...
        'MarkerFaceColor', bearing_colours(bi,:), ...
        'DisplayName', sprintf('%.0f°', batch_test_bearings_deg(bi)));
    hold on;
end
xline(delta_fixed * f_nyq_diag / c_0, '--k', 'HandleVisibility', 'off');
xlabel('\delta / \lambda', 'FontName', fn, 'FontSize', fs_lab);
ylabel('|v_{FD} - v_{true}| / |v_{true}|', 'FontName', fn, 'FontSize', fs_lab);
set(gca, 'FontName', fn, 'FontSize', fs_ax, 'YGrid', 'on', 'Box', 'on');
title('(d) FD velocity error (relative)', 'FontName', fn, 'FontWeight', 'normal');

% (e) Beam pattern: beamwidth + DI on twin axes
subplot(2,3,5);
yyaxis left;
for bi = 1:batch_num_bearings
    bd = bearing_data{bi};
    plot(bd.dl, bd.beamwidth, '-', 'Color', bearing_colours(bi,:), ...
        'LineWidth', lw, 'Marker', bearing_markers{bi}, 'MarkerSize', ms-1, ...
        'MarkerFaceColor', bearing_colours(bi,:), ...
        'DisplayName', sprintf('BW %.0f°', batch_test_bearings_deg(bi)));
    hold on;
end
ylabel('-3 dB beamwidth (deg)', 'FontName', fn, 'FontSize', fs_lab);
yyaxis right;
for bi = 1:batch_num_bearings
    bd = bearing_data{bi};
    plot(bd.dl, bd.di, '--', 'Color', bearing_colours(bi,:), ...
        'LineWidth', lw-0.4, 'Marker', bearing_markers{bi}, 'MarkerSize', ms-3, ...
        'HandleVisibility', 'off');
    hold on;
end
ylabel('DI (dB) - dashed', 'FontName', fn, 'FontSize', fs_lab);
xlabel('\delta / \lambda', 'FontName', fn, 'FontSize', fs_lab);
set(gca, 'FontName', fn, 'FontSize', fs_ax, 'Box', 'on');
title('(e) Beam pattern metrics', 'FontName', fn, 'FontWeight', 'normal');

% (f) Grating lobe count
subplot(2,3,6);
for bi = 1:batch_num_bearings
    bd = bearing_data{bi};
    plot(bd.dl, bd.num_lobes, '-', 'Color', bearing_colours(bi,:), ...
        'LineWidth', lw, 'Marker', bearing_markers{bi}, 'MarkerSize', ms-1, ...
        'MarkerFaceColor', bearing_colours(bi,:), ...
        'DisplayName', sprintf('%.0f°', batch_test_bearings_deg(bi)));
    hold on;
end
xline(delta_fixed * f_nyq_diag / c_0, '--k', 'HandleVisibility', 'off');
xlabel('\delta / \lambda', 'FontName', fn, 'FontSize', fs_lab);
ylabel('Number of lobes', 'FontName', fn, 'FontSize', fs_lab);
ymax_lobes = max(cellfun(@(b) max([b.num_lobes; 1]), bearing_data));
ylim([0.5, ymax_lobes + 0.5]);
set(gca, 'FontName', fn, 'FontSize', fs_ax, 'YGrid', 'on', 'Box', 'on', ...
    'YTick', 1:ymax_lobes);
title('(f) Number of lobes', 'FontName', fn, 'FontWeight', 'normal');

sgtitle(sprintf('Single-AVS Operating Envelope — \\delta = %.3f m, r = %.2f m', ...
    delta_fixed, source_distance_m), 'FontName', fn, 'FontSize', 13);

exportgraphics(gcf, fullfile(results_folder, 'composite_envelope.png'), ...
    'ContentType', 'image');
exportgraphics(gcf, fullfile(results_folder, 'composite_envelope.pdf'), ...
    'ContentType', 'vector');


%% FINAL TEXT SUMMARY %%

fprintf('\n========================================\n');
fprintf('<strong>SWEEP COMPLETE</strong>\n');
fprintf('========================================\n');
fprintf('\nResults: %s\n', results_folder);

% Main results table
fprintf('\n%-9s %-9s %-9s %-13s %-13s %-13s %-7s\n', ...
    'Bearing', 'f (Hz)', 'δ/λ', 'Bracket(r/λ)', 'Radial(m)', 'Angular(°)', 'Sat?');
fprintf('%s\n', repmat('-', 1, 90));
for i = 1:batch_num_tests
    if isnan(batch_res_radial_error(i))
        rad_str = '       NaN';
    else
        rad_str = sprintf('%-13.4f', batch_res_radial_error(i));
    end
    if isnan(batch_res_angular_error(i))
        ang_str = '       NaN';
    else
        ang_str = sprintf('%-13.3f', batch_res_angular_error(i));
    end
    if batch_res_saturated(i)
        sat_str = 'YES';
    else
        sat_str = 'no';
    end
    fprintf('%-9.0f %-9.0f %-9.3f %-13.3f %s %s %-7s\n', ...
        batch_res_bearing_deg(i), batch_res_freq_hz(i), batch_res_dl_ratio(i), ...
        batch_res_bracket_rl(i), rad_str, ang_str, sat_str);
end

% Operational range based on bracket — first δ/λ where saturation kicks in
fprintf('\n<strong>Bracket-defined range-resolution upper edge per bearing:</strong>\n');
for bi = 1:batch_num_bearings
    bd = bearing_data{bi};
    sat_idx = find(bd.saturated, 1, 'first');
    if isempty(sat_idx)
        fprintf('  %.0f°: no saturation in tested range (max δ/λ = %.3f)\n', ...
            batch_test_bearings_deg(bi), bd.dl(end));
    elseif sat_idx == 1
        fprintf('  %.0f°: saturated from lowest tested δ/λ — try lower frequencies\n', ...
            batch_test_bearings_deg(bi));
    else
        fprintf('  %.0f°: saturation onset at δ/λ = %.3f (last clean: %.3f at %.0f Hz)\n', ...
            batch_test_bearings_deg(bi), bd.dl(sat_idx), bd.dl(sat_idx-1), bd.freq(sat_idx-1));
    end
end


%% LOCAL FUNCTIONS %%
% Helpers and MVDR-related functions duplicated from the AVA driver,
% trimmed to single-AVS use. Update both scripts together if any
% shared function changes.

function do_plot_metric(bearing_data, bearings, colours, markers, ...
    field, ylabel_str, use_log, save_path, fn, fs_ax, fs_lab, fs_leg, lw, ms, ...
    delta_fixed, c_0, f_nyq_side, f_nyq_diag)

    figure('Color', 'w', 'Position', [100 100 600 400]);
    for bi = 1:length(bearings)
        bd = bearing_data{bi};
        y_data = bd.(field);
        if use_log
            y_data = max(y_data, eps);
            semilogy(bd.dl, y_data, '-', 'Color', colours(bi,:), ...
                'LineWidth', lw, 'Marker', markers{bi}, 'MarkerSize', ms, ...
                'MarkerFaceColor', colours(bi,:), ...
                'DisplayName', sprintf('%.0f°', bearings(bi)));
        else
            plot(bd.dl, y_data, '-', 'Color', colours(bi,:), ...
                'LineWidth', lw, 'Marker', markers{bi}, 'MarkerSize', ms, ...
                'MarkerFaceColor', colours(bi,:), ...
                'DisplayName', sprintf('%.0f°', bearings(bi)));
        end
        hold on;
    end
    % Reference: diagonal Nyquist line
    xline(delta_fixed * f_nyq_diag / c_0, '--k', 'HandleVisibility', 'off');
    yl = get(gca, 'YLim');
    text(delta_fixed * f_nyq_diag / c_0, yl(2), '  f_{Nyq,diag}', ...
        'FontName', fn, 'FontSize', 9, 'VerticalAlignment', 'top');
    xlabel('\delta / \lambda', 'FontName', fn, 'FontSize', fs_lab);
    ylabel(ylabel_str, 'FontName', fn, 'FontSize', fs_lab);
    legend('Bearing', 'Location', 'best', 'FontName', fn, 'FontSize', fs_leg, 'Box', 'on');
    set(gca, 'FontName', fn, 'FontSize', fs_ax, 'YGrid', 'on', 'Box', 'on');
    exportgraphics(gcf, save_path, 'ContentType', 'image');
end


function plot_2dscan_single(grid_response, x_scan, y_scan, ...
    source_positions, est_pos, vs_centres, bearing_deg, freq)

    figure('Color', 'w', 'Position', [100 100 520 440]);
    imagesc(x_scan, y_scan, grid_response);
    axis xy; axis tight;
    colormap(gca, 'jet');
    cb = colorbar;
    cb.Label.String = 'Normalised response (dB)';
    cb.Label.Interpreter = 'latex';
    xlabel('$x$ (m)', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('$y$ (m)', 'Interpreter', 'latex', 'FontSize', 14);
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 12, ...
        'LineWidth', 1.0, 'Box', 'on', 'Layer', 'top');
    hold on;
    plot(source_positions(:,1), source_positions(:,2), 'kx', ...
        'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', 'True');
    plot(est_pos(1), est_pos(2), 'ko', 'MarkerSize', 10, 'LineWidth', 1.5, ...
        'MarkerFaceColor', 'none', 'DisplayName', 'Estimated');
    ac = mean(vs_centres(1:2,:), 2);
    plot(ac(1), ac(2), 'ks', 'MarkerSize', 9, 'LineWidth', 1, ...
        'MarkerFaceColor', 'w', 'DisplayName', 'Array');
    title(sprintf('Bearing %.0f°, f = %.0f Hz', bearing_deg, freq), ...
        'FontName', 'Times New Roman', 'FontSize', 13);
    legend('Location', 'northeast', 'FontName', 'Times New Roman', ...
        'FontSize', 10, 'Box', 'on');
end


% --- FD VELOCITY ERROR --------------------------------------------------
function [rel_err, ang_err_deg] = compute_velocity_fd_error( ...
    source_pos, mic_positions, vs_centres, freq, c_0, rho_0, N_v)

    % Analytical-field FD velocity error diagnostic.
    %
    % Synthesises the analytical monopole pressure phasor at each mic
    % position, runs the same least-squares finite-difference fit that
    % vs_outputs uses internally, and compares against the analytical
    % Euler-equation velocity at the AVS centre. Both quantities are in
    % the same physical units (no signal-normalisation chain involved),
    % so magnitudes are meaningful and directly comparable to Test B.
    %
    % The diagnostic isolates the geometric question: does the linear-
    % gradient assumption underlying V_3 / V_6 hold across the array?
    % When the relative error is small the assumption holds; when it
    % grows past O(1) the assumption has broken down. Independent of
    % SNR, source amplitude, or signal generation chain.
    %
    % Returns:
    %   rel_err     — |v_fd − v_true| / |v_true| (dimensionless)
    %   ang_err_deg — angle between FD-estimated and true velocity vectors

    omega = 2 * pi * freq;
    k = omega / c_0;

    % Synthesise the analytic monopole pressure at each mic position
    % (no time delay — single-frequency steady-state phasor)
    N_m = N_v * 4;
    p_mic = zeros(N_m, 1);
    for m = 1:N_m
        R = norm(mic_positions(:, m) - source_pos);
        p_mic(m) = exp(-1i * k * R) / R;
    end

    % Per-AVS comparison
    rel_err_per_avs = zeros(N_v, 1);
    ang_err_per_avs = zeros(N_v, 1);
    for vs = 1:N_v
        idx = (vs - 1) * 4 + (1:4);

        % Finite-difference Vx, Vy via the same least-squares fit used
        % in vs_outputs_baseline. Mic positions taken relative to the
        % AVS centre (matches Test B's local-coordinate convention).
        mic_pos_local = (mic_positions(1:2, idx) - vs_centres(1:2, vs)).';
        M = [ones(4, 1), mic_pos_local];
        coeffs = M \ p_mic(idx);
        vx_fd = -coeffs(2) / (1i * omega * rho_0);
        vy_fd = -coeffs(3) / (1i * omega * rho_0);

        % Analytic Vx, Vy at the AVS centre via Euler on the monopole field
        r_vec = vs_centres(:, vs) - source_pos;
        R_centre = norm(r_vec);
        r_hat = r_vec / R_centre;
        common = exp(-1i * k * R_centre) / (1i * omega * rho_0 * R_centre^2);
        vx_true = common * (1 + 1i * k * R_centre) * r_hat(1);
        vy_true = common * (1 + 1i * k * R_centre) * r_hat(2);

        v_fd = [vx_fd; vy_fd];
        v_true = [vx_true; vy_true];

        % Relative error norm
        if norm(v_true) > 0
            rel_err_per_avs(vs) = norm(v_fd - v_true) / norm(v_true);
        else
            rel_err_per_avs(vs) = NaN;
        end

        % Vector angular error — angle between the FD and true velocity
        % vectors. Uses real(v_fd' * v_true) for the cosine since both
        % vectors are complex; phase agreement is implicit in choosing
        % the real part. Clamped for acos numerical safety.
        denom_ang = norm(v_fd) * norm(v_true);
        if denom_ang > 0
            cos_ang = real(v_fd' * v_true) / denom_ang;
            cos_ang = max(-1, min(1, cos_ang));
            ang_err_per_avs(vs) = rad2deg(acos(cos_ang));
        else
            ang_err_per_avs(vs) = NaN;
        end
    end

    rel_err = mean(rel_err_per_avs, 'omitnan');
    ang_err_deg = mean(ang_err_per_avs, 'omitnan');
end


% --- 3 dB RANGE BRACKET ALONG SOURCE-BEARING RAY ------------------------
% Implementation matches Test A's analyse_ray_peak conventions:
%   - Saturation criterion: response within 3 dB of peak at BOTH endpoints
%     (not a fraction-of-search-range heuristic)
%   - When saturated, returned bracket = full search extent (not the
%     contiguous-region width, which can be misleadingly narrow if the
%     above-threshold region is non-contiguous)
%   - Sub-pixel peak refinement via parabolic interpolation
function [bracket_m, bracket_rl, peak_r, ray_response_db, r_vals, saturated] = ...
    compute_range_bracket(r_vs, vs_centres, ray_dir, bin_freqs, c_0, ...
    rho_0, loading, num_bins, N_v, r_min, r_max, lambda, num_r_pts)

    if nargin < 13 || isempty(num_r_pts)
        num_r_pts = 200; %match Test A default
    end

    % Sample the MVDR cost function at num_r_pts points along the ray
    r_vals = linspace(r_min, r_max, num_r_pts);
    response_linear_sum = zeros(num_r_pts, 1);

    for jf = 1:num_bins
        freq = bin_freqs(jf);
        rf = squeeze(r_vs(:,:,jf));
        [ur, er] = eig(rf);
        erv = real(diag(er));
        erv = max(erv, 1e-12);
        lambda_load = loading * max(erv);

        for ir = 1:num_r_pts
            r = r_vals(ir);
            source_pos = [r * ray_dir(1); r * ray_dir(2); 0];
            v_vs = vs_steering_vector_local(vs_centres, source_pos, freq, c_0, rho_0, N_v);
            rxv = (ur * diag(1 ./ (erv + lambda_load)) * ur') * v_vs;
            denominator = v_vs' * rxv;
            if abs(denominator) > 1e-12
                vmvdr = rxv / denominator;
                response_linear_sum(ir) = response_linear_sum(ir) + ...
                    abs(vmvdr' * rf * vmvdr);
            end
        end
    end

    ray_response_db = 10 * log10(response_linear_sum + eps);
    ray_response_db = ray_response_db - max(ray_response_db);

    % Discrete peak
    [pk_val, pk_idx] = max(ray_response_db);

    % Sub-pixel peak refinement via parabolic interpolation (Test A style)
    if pk_idx < 2 || pk_idx > num_r_pts - 1
        peak_r = r_vals(pk_idx);
    else
        dr = r_vals(2) - r_vals(1);
        f_m = ray_response_db(pk_idx - 1);
        f_0 = ray_response_db(pk_idx);
        f_p = ray_response_db(pk_idx + 1);
        denom = f_m - 2 * f_0 + f_p;
        if abs(denom) > eps
            delta_i = (f_m - f_p) / (2 * denom);
        else
            delta_i = 0;
        end
        delta_i = max(-0.5, min(0.5, delta_i));
        peak_r = r_vals(pk_idx) + delta_i * dr;
    end

    % -3 dB range bracket — Test A criterion
    threshold = pk_val - 3;
    above = ray_response_db(:) >= threshold;

    % Saturated when the response is within 3 dB of the peak at BOTH the
    % first and last sample of the ray. This is the operationally
    % meaningful "cost function spilled out of the search range" check.
    saturated = above(1) && above(end);

    if saturated
        bracket_m = r_vals(end) - r_vals(1);
    else
        d_above = diff([0; above; 0]);
        starts = find(d_above == 1);
        ends_arr = find(d_above == -1) - 1;
        region = find(starts <= pk_idx & ends_arr >= pk_idx, 1);
        if isempty(region)
            bracket_m = NaN;
        else
            bracket_m = r_vals(ends_arr(region)) - r_vals(starts(region));
        end
    end

    bracket_rl = bracket_m / lambda;
end


% --- TERNARY HELPER -----------------------------------------------------
function s = ternary_str(cond, str_true, str_false)
    if cond
        s = str_true;
    else
        s = str_false;
    end
end


% --- MVDR / SIGNAL PROCESSING (single-AVS subset) -----------------------

function snapshots = make_snapshots(tx_vs, size_fft, overlap, window)
    step = round(overlap * size_fft);
    num_snap = floor((size(tx_vs, 2) - size_fft) / step) + 1;
    start_idx = 1 + (0:(num_snap - 1)) * step;
    idx_matrix = start_idx + (0:size_fft - 1)';
    snapshots = tx_vs(:, idx_matrix);
    snapshots = reshape(snapshots, size(tx_vs, 1), size_fft, num_snap);
    snapshots = snapshots .* reshape(window, [1, size_fft, 1]);
end

function r_vs = create_vs_csm(snapshots_vs, bin_index, delta, bin_freqs, ...
    rho_0, N_v, num_bins, c_0)

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
            [p_vs, vx_vs, vy_vs] = vs_outputs_local(fx_snap, delta, freq, rho_0, N_v, c_0);
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
        lambda_load = loading * max(erv);

        for n = 1:num_points
            source_pos = candidate_points(n,:).';
            v_vs = vs_steering_vector_local(vs_centres, source_pos, freq, c_0, rho_0, N_v);
            rxv = (ur * diag(1 ./ (erv + lambda_load)) * ur') * v_vs;
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

function [p_vs, vx_vs, vy_vs] = vs_outputs_local(tx_freq, delta, freq, ...
    rho_0, N_v, c_0)
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

function v_vs = vs_steering_vector_local(vs_centres, source_pos, freq, ...
    c_0, rho_0, N_v)
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

function [beam_metrics, fig_handle] = compute_nearfield_beam_pattern(r_vs, ...
    vs_centres, array_centre, source_positions, bin_freqs, c_0, rho_0, ...
    loading, num_bins, N_v, radius, plot_title, grating_threshold_db)

    % grating_threshold_db is the dB drop from the mainlobe peak below
    % which a competing peak still counts as a grating lobe. Default 3 dB
    % matches the standard "another peak as strong as the mainlobe"
    % criterion. Pass 6 for a more lenient counter.
    if nargin < 13 || isempty(grating_threshold_db)
        grating_threshold_db = 3.0;
    end

    num_angles = 360;
    theta_deg = linspace(0, 360, num_angles + 1);
    theta_deg = theta_deg(1:end-1);
    theta_rad = deg2rad(theta_deg);
    mvdr_responses = zeros(num_angles, num_bins);

    for jf = 1:num_bins
        freq = bin_freqs(jf);
        rf = squeeze(r_vs(:,:,jf));
        [ur, er] = eig(rf);
        erv = real(diag(er));
        erv = max(erv, 1e-12);
        lambda_load = loading * max(erv);

        for ai = 1:num_angles
            x_pos = array_centre(1) + radius * sin(theta_rad(ai));
            y_pos = array_centre(2) + radius * cos(theta_rad(ai));
            source_pos = [x_pos; y_pos; 0];
            v_vs = vs_steering_vector_local(vs_centres, source_pos, freq, c_0, rho_0, N_v);
            rxv = (ur * diag(1 ./ (erv + lambda_load)) * ur') * v_vs;
            denominator = v_vs' * rxv;
            if abs(denominator) > 1e-12
                vmvdr = rxv / denominator;
                mvdr_responses(ai, jf) = abs(vmvdr' * rf * vmvdr);
            end
        end
    end

    beam_linear = sum(mvdr_responses, 2);
    beam_pattern = 10 * log10(beam_linear + eps);
    beam_pattern = beam_pattern - max(beam_pattern);
    bp_col = beam_pattern(:);

    % --- Beamwidth, mainlobe region, sidelobe ---
    beam_metrics = struct();
    [~, peak_idx] = max(beam_pattern);
    threshold_3db = max(beam_pattern) - 3;
    above = bp_col >= threshold_3db;
    d_above = diff([0; above; 0]);
    starts = find(d_above == 1);
    ends_arr = find(d_above == -1) - 1;

    beam_metrics.beamwidth_3dB = NaN;
    main_lobe_region = find(starts <= peak_idx & ends_arr >= peak_idx, 1);
    if ~isempty(main_lobe_region)
        bw = theta_deg(ends_arr(main_lobe_region)) - theta_deg(starts(main_lobe_region));
        if bw < 0, bw = bw + 360; end
        beam_metrics.beamwidth_3dB = bw;
    end

    beam_metrics.sidelobe_level = NaN;
    if any(~above)
        beam_metrics.sidelobe_level = max(bp_col(~above));
    end

    beam_linear_full = 10.^(bp_col / 10);
    integral_val = trapz(theta_rad(:), beam_linear_full);
    beam_metrics.directivity_index = 10 * log10(2 * pi / integral_val);


    % --- GRATING LOBE DETECTION ---
    % Three-step approach to avoid the mainlobe-as-grating-lobe artefact
    % that the previous (single-pass local-max) detector suffered from:
    %
    %   (1) Find the mainlobe peak via global max (with parabolic refinement
    %       skipped — angular peak resolution is 1° already)
    %   (2) Mask out the mainlobe lobe region — set the response to -Inf
    %       within ±half_beamwidth of the mainlobe angle. This guarantees
    %       no peak detected later can be the mainlobe itself
    %   (3) Find the strongest remaining peak above the threshold; if any,
    %       record it as the grating lobe
    %
    % This is more robust than the previous "find all local maxima then
    % filter by separation > half-bw" approach, which was failing when
    % the cardioid pattern degraded into a multi-lobe pattern with several
    % peaks at near-zero dB (the back lobe and other side lobes were all
    % being flagged as grating lobes of the mainlobe).

    % Mainlobe = global max
    [mainlobe_level, mainlobe_idx] = max(bp_col);
    mainlobe_angle_deg = theta_deg(mainlobe_idx);

    if isnan(beam_metrics.beamwidth_3dB)
        half_bw_deg = 30; %conservative fallback if mainlobe couldn't be measured
    else
        half_bw_deg = beam_metrics.beamwidth_3dB / 2;
    end

    % Mask out the mainlobe lobe region — angular distance from mainlobe
    % within half-beamwidth gets set to -Inf so it cannot be re-detected.
    dtheta_to_main = abs(theta_deg - mainlobe_angle_deg);
    dtheta_to_main = min(dtheta_to_main, 360 - dtheta_to_main); %wrap to [0, 180]
    masked_pattern = bp_col;
    masked_pattern(dtheta_to_main(:) <= half_bw_deg) = -Inf;

    % Search for the strongest remaining peak. We look for a local maximum
    % on the masked pattern; if there is a peak above the threshold, it's
    % a grating lobe.
    finite_mask = isfinite(masked_pattern);
    if any(finite_mask)
        [grating_level, grating_idx] = max(masked_pattern);
        grating_angle_deg = theta_deg(grating_idx);
        grating_sep_deg = abs(grating_angle_deg - mainlobe_angle_deg);
        grating_sep_deg = min(grating_sep_deg, 360 - grating_sep_deg);

        % Apply the threshold criterion — must be within `grating_threshold_db`
        % of the mainlobe (which sits at 0 dB after normalisation)
        is_grating = (grating_level >= -grating_threshold_db);
    else
        is_grating = false;
        grating_level = NaN;
        grating_angle_deg = NaN;
        grating_sep_deg = NaN;
    end

    if is_grating
        beam_metrics.num_lobes = 2;
        beam_metrics.grating_lobe_present = true;
        beam_metrics.grating_lobe_angle_deg = grating_angle_deg;
        beam_metrics.grating_lobe_separation_deg = grating_sep_deg;
        beam_metrics.grating_lobe_level_db = grating_level;
    else
        beam_metrics.num_lobes = 1;
        beam_metrics.grating_lobe_present = false;
        beam_metrics.grating_lobe_angle_deg = NaN;
        beam_metrics.grating_lobe_separation_deg = NaN;
        beam_metrics.grating_lobe_level_db = NaN;
    end


    % --- PLOT ---
    fig_handle = figure('Name', plot_title, 'Position', [100, 100, 600, 600]);
    polarplot(theta_rad, beam_pattern, 'b-', 'LineWidth', 2);
    ax = gca;
    ax.ThetaDir = 'clockwise';
    ax.ThetaZeroLocation = 'top';
    rlim([min(beam_pattern), 0]);

    hold on;

    % Mark the grating lobe if present
    if is_grating
        polarplot(theta_rad(grating_idx), bp_col(grating_idx), 'rs', ...
            'MarkerSize', 8, 'LineWidth', 1.5, 'MarkerFaceColor', 'none');
    end

    % Source line
    if ~isempty(source_positions)
        for src = 1:size(source_positions, 1)
            sv = source_positions(src, :)' - array_centre;
            sa = atan2(sv(1), sv(2));
            polarplot([sa, sa], rlim, '--r', 'LineWidth', 1.5);
        end
    end
    title(plot_title, 'FontName', 'Times New Roman', 'FontSize', 13);

    fprintf('  Beamwidth: %.1f° | Sidelobe: %.2f dB | DI: %.2f dB\n', ...
        beam_metrics.beamwidth_3dB, beam_metrics.sidelobe_level, ...
        beam_metrics.directivity_index);
    if beam_metrics.grating_lobe_present
        fprintf('  GRATING LOBE: %d total lobes, strongest at %.1f° (Δ%.1f°), %.2f dB\n', ...
            beam_metrics.num_lobes, ...
            beam_metrics.grating_lobe_angle_deg, ...
            beam_metrics.grating_lobe_separation_deg, ...
            beam_metrics.grating_lobe_level_db);
    else
        fprintf('  No grating lobe (within %.1f dB of mainlobe)\n', grating_threshold_db);
    end
end


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
    x_pts = round(2 * max_extent / grid_resolution) + 1;
    y_pts = round(2 * max_extent / grid_resolution) + 1;
    x_scan = linspace(-max_extent, max_extent, x_pts);
    y_scan = linspace(array_centre_y - max_extent, array_centre_y + max_extent, y_pts);
    [X_grid, Y_grid] = meshgrid(x_scan, y_scan);
    candidate_points = [X_grid(:), Y_grid(:), zeros(numel(X_grid), 1)];
    grid_res = mean([x_scan(2)-x_scan(1), y_scan(2)-y_scan(1)]);
end


function [est_x, est_y, refined_db] = refine_peak_2d(grid_response, x_scan, y_scan)
    [~, max_idx] = max(grid_response(:));
    [iy_pk, ix_pk] = ind2sub(size(grid_response), max_idx);
    est_x_d = x_scan(ix_pk); est_y_d = y_scan(iy_pk);
    dx = x_scan(2) - x_scan(1); dy = y_scan(2) - y_scan(1);

    if ix_pk < 2 || ix_pk > length(x_scan) - 1 || ...
       iy_pk < 2 || iy_pk > length(y_scan) - 1
        est_x = est_x_d; est_y = est_y_d;
        refined_db = grid_response(iy_pk, ix_pk);
        return;
    end

    fx_m = grid_response(iy_pk, ix_pk - 1);
    fx_0 = grid_response(iy_pk, ix_pk);
    fx_p = grid_response(iy_pk, ix_pk + 1);
    denom_x = fx_m - 2 * fx_0 + fx_p;
    if abs(denom_x) > eps, dix = (fx_m - fx_p)/(2*denom_x); else, dix = 0; end

    fy_m = grid_response(iy_pk - 1, ix_pk);
    fy_p = grid_response(iy_pk + 1, ix_pk);
    denom_y = fy_m - 2 * fx_0 + fy_p;
    if abs(denom_y) > eps, diy = (fy_m - fy_p)/(2*denom_y); else, diy = 0; end

    dix = max(-0.5, min(0.5, dix));
    diy = max(-0.5, min(0.5, diy));
    est_x = est_x_d + dix * dx;
    est_y = est_y_d + diy * dy;
    refined_db = fx_0 - (fx_m - fx_p)^2 / (8 * denom_x);
end