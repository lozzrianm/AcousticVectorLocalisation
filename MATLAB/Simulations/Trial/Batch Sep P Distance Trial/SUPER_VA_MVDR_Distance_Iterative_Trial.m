clc; close all;

% L Marshall - CORRECTED Iterative Range-Adaptive MVDR Beamforming
% 
% KEY INSIGHT: The issue with the original approach was applying the correction
% in the wrong direction. For spherical waves at short ranges:
% - Velocity magnitude is LARGER than plane wave prediction
% - The correction factor should be (1 + 1/(ikr)), not (ikr)/(1 + ikr)
%
% PHYSICS: For a monopole at range r:
%   v_spherical = v_plane × (1 + 1/(ikr))
% where:
%   - kr << 1 (near field): factor ≈ 1/(ikr) >> 1 (velocity enhanced)
%   - kr >> 1 (far field): factor ≈ 1 (plane wave limit)

%% TRIAL PARAMETERS %%

max_iterations = 3;
convergence_threshold = 0.001;  % Tighter threshold (1mm)

%% CHECK FOR BATCH MODE %%

batch_mode = exist('batch_csv_input', 'var');

if batch_mode
    csv_filename = batch_csv_input;
    if ~exist('batch_save_figures', 'var')
        set(0, 'DefaultFigureVisible', 'off');
    end
else
    csv_filename = 'generatedsignal_avs.csv';
end

%% DEFINE INPUT VARIABLES %%

c_0 = 340;
rho_0 = 1.02;

if exist('batch_test_source_x', 'var') && exist('batch_test_source_y', 'var')
    source_positions = [batch_test_source_x, batch_test_source_y, 0];
else
    source_positions = [-0.7071, 0.7071, 0];
end

if exist('batch_test_freq', 'var')
    source_frequencies = batch_test_freq;
else
    source_frequencies = [1500];
end

source_amplitudes = [0.01];
num_sources = size(source_positions, 1);

N_a = 1;
N_v = 1;
d_y = 0.1;

if exist('batch_test_delta', 'var')
    delta = batch_test_delta;
else
    delta = 0.042;
end

x_a1 = 0; z_a1 = 0; y_a1 = 0;
A1_coord = [x_a1; y_a1; z_a1];

lambda = c_0 / source_frequencies(1, :);

freq_limit = 3000;
overlap = 0.5;
loading = 1e-4;

if exist('batch_x_scan_points', 'var')
    x_scan_points = batch_x_scan_points;
    y_scan_points = batch_y_scan_points;
    x_margin = batch_x_margin;
    y_margin = batch_y_margin;
else
    x_scan_points = 200;
    y_scan_points = 200;
    x_margin = 3;
    y_margin = 2;
end

fprintf('\n<strong>CORRECTED ITERATIVE RANGE-ADAPTIVE BEAMFORMING</strong>\n');
fprintf('  Source: (%.3f, %.3f) m\n', source_positions(1, 1:2));
fprintf('  Frequency: %.0f Hz (lambda = %.4f m)\n', source_frequencies(1), lambda);
fprintf('  True range: %.4f m (%.3f lambda)\n', ...
    norm(source_positions(1, 1:2)), norm(source_positions(1, 1:2))/lambda);

%% LOAD AND PROCESS SIGNAL %%

data_table = readtable(csv_filename);
data = table2array(data_table);
time = data(:, 1);
N = size(data, 1);
N_ma = N_v * 4;
N_m = N_a * N_ma;

tx_vs_arrays = cell(N_a, 1);
for array = 1:N_a
    tx_vs = zeros(N_ma, N);
    for idx_local = 1:N_ma
        monopole_col = 1 + (array-1)*N_ma + idx_local;
        sinusoidal_col = 1 + N_m + (array-1)*N_ma + idx_local;
        tx_vs(idx_local, :) = data(:, monopole_col) + data(:, sinusoidal_col);
    end
    tx_vs_arrays{array} = tx_vs;
end

%% ARRAY GEOMETRY %%

array_centres = A1_coord;
if N_v > 1
    array_centres(2, :) = array_centres(2, :) - (N_v-1)*d_y/2;
end

vs_centres_arrays = cell(N_a, 1);
mic_offsets = delta/2 * [-1, 1, 1, -1; -1, -1, 1, 1; 0, 0, 0, 0];

for array = 1:N_a
    vs_centres = zeros(3, N_v);
    for vs = 1:N_v
        vs_centres(:, vs) = array_centres(:, array) + [0; (vs-1)*d_y; 0];
    end
    vs_centres_arrays{array} = vs_centres;
end

%% CONVERT TO FREQUENCY DOMAIN %%

d_t = time(2) - time(1);
F_s = 1/d_t;
f = F_s * (0:floor(N/2)) / N;

if batch_mode
    freq_limit = max(freq_limit, source_frequencies(1) * 1.5);
end
idx_limit = f <= freq_limit;

tx_freq_arrays = cell(N_a, 1);
for array = 1:N_a
    tx_freq_arrays{array} = zeros(N_ma, sum(idx_limit));
    for sensor = 1:N_ma
        Y = fft(tx_vs_arrays{array}(sensor, :));
        Y = Y(1:floor(N/2)+1);
        tx_freq_arrays{array}(sensor, :) = Y(idx_limit);
    end
end

%% DEFINE GRID %%

x_scan = linspace(x_a1 - x_margin, x_a1 + x_margin, x_scan_points);
y_scan = linspace(y_a1 - y_margin, y_a1 + y_margin, y_scan_points);
[X_grid, Y_grid] = meshgrid(x_scan, y_scan);
Z_grid = zeros(size(X_grid));
candidate_points = [X_grid(:), Y_grid(:), Z_grid(:)];
grid_res = mean([x_scan(2)-x_scan(1), y_scan(2)-y_scan(1)]);

%% CREATE FREQUENCY BINS %%

max_binwidth = 1/(8*d_y*(N_ma-1)/c_0);
size_fft = floor(F_s/max_binwidth);
fft_vec = F_s * (0:(size_fft-1)) / size_fft;
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
window = hanning(size_fft)';

%% CREATE SNAPSHOTS %%

snapshots_vs_arrays = cell(N_a, 1);
for array = 1:N_a
    snapshots_vs_arrays{array} = make_snapshots(tx_vs_arrays{array}, size_fft, overlap, window);
end

%% ITERATIVE BEAMFORMING %%

fprintf('\n========================================\n');
fprintf('<strong>STARTING ITERATIONS</strong>\n');
fprintf('========================================\n');

iteration_history = struct('est_positions', cell(max_iterations + 1, 1), ...
                          'est_ranges', cell(max_iterations + 1, 1), ...
                          'response_db', cell(max_iterations + 1, 1), ...
                          'radial_errors', zeros(max_iterations + 1, 1));

for iter = 0:max_iterations
    
    if iter == 0
        fprintf('\n<strong>ITERATION 0: Baseline</strong>\n');
        estimated_ranges = [];
    else
        fprintf('\n<strong>ITERATION %d</strong>\n', iter);
        estimated_ranges = iteration_history(iter).est_ranges;
        k = 2*pi*source_frequencies(1)/c_0;
        kr = k * estimated_ranges;
        correction_mag = abs(1 + 1/(1i*k*estimated_ranges));
        fprintf('  Using r_est = %.4f m (%.3f lambda, kr = %.2f)\n', ...
            estimated_ranges, estimated_ranges/lambda, kr);
        fprintf('  Correction magnitude: %.3f (%.0f%% velocity boost)\n', ...
            correction_mag, 100*(correction_mag - 1));
    end
    
    %% CREATE CSM WITH CORRECTED VELOCITY ESTIMATES %%
    
    r_vs_arrays = cell(N_a, 1);
    
    for array = 1:N_a
        r_vs_arrays{array} = create_vs_csm_corrected(snapshots_vs_arrays{array}, ...
            bin_index, delta, bin_freqs, rho_0, N_v, num_bins, c_0, estimated_ranges);
    end
    
    %% PERFORM MVDR BEAMFORMING %%
    
    response_db_current = mvdr_vs_beamforming(r_vs_arrays{1}, vs_centres_arrays{1}, ...
        candidate_points, bin_freqs, c_0, rho_0, loading, num_bins, N_v);
    
    %% LOCALISE SOURCE %%
    
    grid_response = reshape(response_db_current, length(y_scan), length(x_scan));
    [~, max_idx] = max(grid_response(:));
    [y_idx, x_idx] = ind2sub(size(grid_response), max_idx);
    est_position = [X_grid(y_idx, x_idx), Y_grid(y_idx, x_idx)];
    
    array_centre_2d = vs_centres_arrays{1}(1:2, 1);
    est_range = norm(est_position - array_centre_2d');
    
    true_pos = source_positions(1, 1:2);
    radial_error = norm(true_pos - est_position);
    
    iteration_history(iter + 1).est_positions = est_position;
    iteration_history(iter + 1).est_ranges = est_range;
    iteration_history(iter + 1).response_db = response_db_current;
    iteration_history(iter + 1).radial_errors = radial_error;
    
    fprintf('  Estimated: (%.4f, %.4f) m, range: %.4f m\n', est_position, est_range);
    fprintf('  Error: %.4f m\n', radial_error);
    
    if iter > 0
        pos_change = norm(est_position - iteration_history(iter).est_positions);
        error_change = iteration_history(iter).radial_errors - radial_error;
        fprintf('  Position change: %.5f m\n', pos_change);
        fprintf('  Error improvement: %.5f m (%.1f%%)\n', error_change, ...
            100*error_change/iteration_history(iter).radial_errors);
        
        if pos_change < convergence_threshold
            fprintf('  <strong>Converged!</strong>\n');
            final_iteration = iter;
            break;
        end
    end
    
    final_iteration = iter;
end

fprintf('\n========================================\n');
fprintf('<strong>COMPLETE</strong>\n');
fprintf('========================================\n');

%% RESULTS SUMMARY %%

fprintf('\n%-10s %-20s %-15s\n', 'Iteration', 'Position (m)', 'Error (m)');
fprintf('%s\n', repmat('-', 1, 50));

for iter = 0:final_iteration
    pos = iteration_history(iter + 1).est_positions;
    error = iteration_history(iter + 1).radial_errors;
    fprintf('%-10d (%.4f, %.4f)    %.4f\n', iter, pos, error);
end

baseline_error = iteration_history(1).radial_errors;
final_error = iteration_history(final_iteration + 1).radial_errors;
improvement = baseline_error - final_error;

fprintf('\nBaseline error: %.4f m\n', baseline_error);
fprintf('Final error: %.4f m\n', final_error);
fprintf('Improvement: %.4f m (%.1f%%)\n', improvement, 100*improvement/baseline_error);

%% STORE FOR BATCH COMPATIBILITY %%

all_results = struct('est_positions', iteration_history(final_iteration + 1).est_positions, ...
                    'true_positions', source_positions(1, 1:2), ...
                    'radial_error', final_error, ...
                    'angular_error', NaN, ...
                    'peak_response', max(iteration_history(final_iteration + 1).response_db));

true_angle = atan2(source_positions(1, 2), source_positions(1, 1));
est_angle = atan2(iteration_history(final_iteration + 1).est_positions(2), ...
                 iteration_history(final_iteration + 1).est_positions(1));
angular_error_deg = rad2deg(abs(true_angle - est_angle));
if angular_error_deg > 180
    angular_error_deg = 360 - angular_error_deg;
end
all_results.angular_error = angular_error_deg;

beam_metrics = struct('beamwidth_3dB', NaN, 'sidelobe_level', NaN, 'directivity_index', NaN);

if batch_mode
    set(0, 'DefaultFigureVisible', 'on');
end


%% FUNCTION DEFINITIONS %%

function snapshots = make_snapshots(tx_vs, size_fft, overlap, window)
    step = round(overlap * size_fft);
    num_snap = floor((size(tx_vs, 2) - size_fft) / step) + 1;
    start_idx = 1 + (0:(num_snap-1)) * step;
    idx_matrix = start_idx + (0:size_fft-1)';
    snapshots = tx_vs(:, idx_matrix);
    snapshots = reshape(snapshots, size(tx_vs, 1), size_fft, num_snap);
    snapshots = snapshots .* reshape(window, [1, size_fft, 1]);
end


function r_vs = create_vs_csm_corrected(snapshots_vs, bin_index, delta, bin_freqs, ...
    rho_0, N_v, num_bins, c_0, estimated_ranges)
    
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
            
            [p_vs, vx_vs, vy_vs] = vs_outputs_corrected(fx_snap, delta, freq, ...
                rho_0, N_v, c_0, estimated_ranges);
            
            rho_c = rho_0 * c_0;
            vx_vs = rho_c * vx_vs;
            vy_vs = rho_c * vy_vs;
            
            vs_signal = zeros(signal_dim, 1);
            for n = 1:N_v
                idx = (n-1)*3 + (1:3);
                vs_signal(idx) = [p_vs(n); vx_vs(n); vy_vs(n)];
            end
            
            R_freq = R_freq + (vs_signal * vs_signal');
        end
        
        r_vs(:, :, jf) = R_freq / num_snap;
    end
end


function [p_vs, vx_vs, vy_vs] = vs_outputs_corrected(tx_freq, delta, freq, rho_0, ...
    N_v, c_0, estimated_ranges)
    
    omega = 2*pi*freq;
    k = omega / c_0;
    p_vs = zeros(N_v, 1);
    vx_vs = zeros(N_v, 1);
    vy_vs = zeros(N_v, 1);
    
    for vs = 1:N_v
        idx = (vs-1)*4 + (1:4);
        p0 = tx_freq(idx(1));
        p1 = tx_freq(idx(2));
        p2 = tx_freq(idx(3));
        p3 = tx_freq(idx(4));
        
        mic_positions = [
            -delta/2, -delta/2;
             delta/2, -delta/2;
             delta/2,  delta/2;
            -delta/2,  delta/2
        ];
        
        M = [ones(4, 1), mic_positions];
        p_vector = [p0; p1; p2; p3];
        
        coeffs = M \ p_vector;
        p_vs(vs) = coeffs(1);
        dpdx = coeffs(2);
        dpdy = coeffs(3);
        
        % Standard plane wave velocity estimate
        vx_plane = -dpdx / (1i * omega * rho_0);
        vy_plane = -dpdy / (1i * omega * rho_0);
        
        if ~isempty(estimated_ranges) && estimated_ranges > 0
            % CORRECTED spherical wave factor
            % For spherical waves: v = v_plane × (1 + 1/(ikr))
            % This INCREASES velocity at short ranges (near field)
            correction = 1 + 1/(1i*k*estimated_ranges);
            
            vx_vs(vs) = vx_plane * correction;
            vy_vs(vs) = vy_plane * correction;
        else
            vx_vs(vs) = vx_plane;
            vy_vs(vs) = vy_plane;
        end
    end
end


function [response_db, mean_condition_number] = mvdr_vs_beamforming(r_vs, vs_centres, ...
    candidate_points, bin_freqs, c_0, rho_0, loading, num_bins, N_v)
    
    num_points = size(candidate_points, 1);
    mvdr_responses = zeros(num_points, num_bins);
    
    for jf = 1:num_bins 
        freq = bin_freqs(jf);
        rf = squeeze(r_vs(:, :, jf));
        
        [ur, er] = eig(rf);
        erv = real(diag(er));
        erv = max(erv, 1e-12);

        for n = 1:num_points
            source_pos = candidate_points(n, :).';
            v_vs = vs_steering_vector(vs_centres, source_pos, freq, c_0, rho_0, N_v);
            
            lambda = loading * max(erv);
            rxv = (ur * diag(1./(erv + lambda)) * ur') * v_vs;
            denominator = v_vs' * rxv;
            
            if abs(denominator) > 1e-12
                vmvdr = rxv / denominator;
                mvdr_responses(n, jf) = abs(vmvdr' * rf * vmvdr);
            end
        end
    end
    
    responses_sum = sum(mvdr_responses, 2);
    response_db = 10*log10(responses_sum + eps);
    response_db = response_db - max(response_db);
    mean_condition_number = NaN;
end


function v_vs = vs_steering_vector(vs_centres, source_pos, freq, c_0, rho_0, N_v)
    k = 2*pi*freq/c_0;
    omega = 2*pi*freq;
    
    v_vs = zeros(3*N_v, 1);
    
    for n = 1:N_v
        r_vec = vs_centres(:, n) - source_pos;
        r = norm(r_vec);
        r_hat = r_vec / r;
        
        p_steer = exp(-1i*k*r) / r;
        
        common_term = exp(-1i*k*r) / (1i*omega*rho_0*r^2);
        velocity_factor = (1 + 1i*k*r);
        
        vx_steer = common_term * velocity_factor * r_hat(1);
        vy_steer = common_term * velocity_factor * r_hat(2);
        
        rho_c = rho_0*c_0;
        vx_steer = rho_c * vx_steer;
        vy_steer = rho_c * vy_steer;
        
        idx = (n-1)*3 + (1:3);
        v_vs(idx) = [p_steer; vx_steer; vy_steer];
    end
end