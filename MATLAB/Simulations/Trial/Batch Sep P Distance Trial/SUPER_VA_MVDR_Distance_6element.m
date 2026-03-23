clc; close all; %clear all;

% L Marshall 27/11/2025
% Variable Number of Sources Localisation via an Acoustic Vector Array
% MVDR method - 6-ELEMENT VERSION with Distance Sweep support

%% CHECK FOR BATCH MODE PARAMETERS %%

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

% Physical Environment Parameters
c_0 = 340;
rho_0 = 1.02;

% Source position - from batch script if available
if exist('batch_test_source_x', 'var') && exist('batch_test_source_y', 'var')
    source_positions = [batch_test_source_x, batch_test_source_y, 0];
else
    source_positions = [-0.7071, 0.7071, 0];
end

% Source frequencies
if exist('batch_test_freq', 'var')
    source_frequencies = batch_test_freq;
else
    source_frequencies = [1500];
end

source_amplitudes = [0.01];

num_sources = size(source_positions, 1);

% Vector Sensor Array Parameters
N_a = 1;
N_v = 1;
d_y = 0.1;

if exist('batch_test_delta', 'var')
    delta = batch_test_delta;
else
    delta = 0.042;
end

% Array centre
x_a1 = 0;
z_a1 = 0;
y_a1 = 0;
A1_coord = [x_a1;y_a1;z_a1];

lambda = c_0/source_frequencies(1,:);
a_spacing = sqrt(((10*lambda)^2)/2);
x_a2 = x_a1 - a_spacing;
z_a2 = 0;
y_a2 = y_a1 - a_spacing;
A2_coord = [x_a2; y_a2; z_a2];

% Sampling Characteristics
freq_limit = 3000;
overlap = 0.5;
loading = 1e-2;

% Grid Search Parameters - can be overridden by batch script
if exist('batch_x_scan_points', 'var')
    x_scan_points = batch_x_scan_points;
else
    x_scan_points = 200;
end

if exist('batch_y_scan_points', 'var')
    y_scan_points = batch_y_scan_points;
else
    y_scan_points = 200;
end

if exist('batch_x_margin', 'var')
    x_margin = batch_x_margin;
else
    x_margin = 3;
end

if exist('batch_y_margin', 'var')
    y_margin = batch_y_margin;
else
    y_margin = 2;
end

if batch_mode
    fprintf('\n<strong>BATCH MODE ACTIVE - 6-ELEMENT Distance Sweep</strong>\n');
    fprintf('  CSV input: %s\n', csv_filename);
    fprintf('  Source: (%.3f, %.3f) m\n', source_positions(1,1:2));
    fprintf('  Frequency: %.0f Hz\n', source_frequencies(1));
    fprintf('  Delta: %.4f m\n', delta);
    fprintf('  Grid: %d x %d points, margins: %.3f x %.3f m\n', ...
        x_scan_points, y_scan_points, x_margin, y_margin);
end

fprintf('\n<strong>6-ELEMENT STEERING VECTOR</strong>\n');
fprintf('Elements: [p0, p1, p2, p3, vx, vy] per vector sensor\n\n');

%% LOAD GENERATED SIGNAL FROM CSV %%

fprintf('Loading signal data from: %s\n', csv_filename);

data_table = readtable(csv_filename);
data = table2array(data_table);
time = data(:,1);
[numRows, numCols] = size(data);

N = numRows;
N_ma = N_v * 4;
N_m = N_a * N_ma;

fprintf('Loaded %d samples from %d microphones (%d arrays, %d mics each)\n', ...
    N, N_m, N_a, N_ma);

%% COMBINE SIGNAL COMPONENTS FOR PROCESSING %%

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

%% ARRAY GEOMETRY SET UP %%

array_centres = [A1_coord, A2_coord];

if N_v > 1
    array_centres(2, :) = array_centres(2, :) - (N_v-1)*d_y/2;
end

vs_centres_arrays = cell(N_a, 1);
mic_positions = zeros(3, N_m);
mic_offsets = delta/2 * [-1, 1, 1, -1; -1, -1, 1, 1; 0, 0, 0, 0];

% Store individual mic positions for each VS
mic_positions_per_vs = cell(N_a, 1);

for array = 1:N_a
    vs_centres = zeros(3, N_v);
    mic_pos_array = zeros(3, 4, N_v);
    
    for vs = 1:N_v
        vs_centres(:, vs) = array_centres(:, array) + [0; (vs-1)*d_y; 0];
        
        idx = (array-1)*N_ma + (vs-1)*4 + (1:4);
        mic_positions(:, idx) = vs_centres(:, vs) + mic_offsets;
        
        mic_pos_array(:, :, vs) = vs_centres(:, vs) + mic_offsets;
    end
    
    vs_centres_arrays{array} = vs_centres;
    mic_positions_per_vs{array} = mic_pos_array;
end

if N_a > 1
    fprintf('Array separation: %.3f m\n', norm(diff(array_centres, 1, 2)));
end

%% CONVERT TO FREQUENCY DOMAIN %%

fprintf('Converting to frequency domain...\n');

d_t = time(2) - time(1);
F_s = 1/d_t;
f = F_s * (0:floor(N/2)) / N;

if batch_mode
    freq_limit = max(freq_limit, source_frequencies(1) * 1.5);
end
idx_limit = f <= freq_limit;

tx_vs_arrays_time = tx_vs_arrays;
tx_freq_arrays = cell(N_a, 1);

for array = 1:N_a
    tx_freq_arrays{array} = zeros(N_ma, sum(idx_limit));
    
    for sensor = 1:N_ma
        Y = fft(tx_vs_arrays_time{array}(sensor, :));
        Y = Y(1:floor(N/2)+1);
        tx_freq_arrays{array}(sensor, :) = Y(idx_limit);
    end
end

%% DEFINING THE GRID SEARCH AREA %%

fprintf('Defining search area...\n');

x_arrays_mean = mean(array_centres(1, :));
y_arrays_mean = mean(array_centres(2, :));

if N_a == 1
    x_scan = linspace(x_a1 - x_margin, x_a1 + x_margin, x_scan_points);
    y_scan = linspace(y_a1 - y_margin, y_a1 + y_margin, y_scan_points);
    z_scan = 0;
else
    x_scan = linspace(x_arrays_mean - x_margin, x_arrays_mean + x_margin, x_scan_points);
    y_scan = linspace(y_arrays_mean - y_margin, y_arrays_mean + y_margin, y_scan_points);
    z_scan = 0;
end

[X_grid, Y_grid] = meshgrid(x_scan, y_scan);
Z_grid = zeros(size(X_grid));

candidate_points = [X_grid(:), Y_grid(:), Z_grid(:)];

x_res = x_scan(2) - x_scan(1);
y_res = y_scan(2) - y_scan(1);
grid_res = mean([x_res, y_res]);

fprintf('Search grid: %d x %d points, resolution: %.4f m\n', ...
    length(x_scan), length(y_scan), grid_res);

%% CREATING FREQUENCY BINS FOR BROADBAND ANALYSIS %%

fprintf('Creating frequency bins...\n');

max_binwidth = 1/(8*d_y*(N_ma-1)/c_0);
size_fft = floor(F_s/max_binwidth);

delta_f = F_s / N;
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

fprintf('Using %d frequency bins (%.1f - %.1f Hz)\n', ...
    num_bins, min(bin_freqs), max(bin_freqs));

num_snaps = 8*N_m;
num_samps = (num_snaps+1)*size_fft*overlap;
window = hanning(size_fft)';

%% CREATE AVS SNAPSHOTS AND MAKE CSM %%

fprintf('Creating snapshots and cross-spectral matrix (6-element)...\n');

snapshots_vs_arrays = cell(N_a, 1);
r_vs_arrays = cell(N_a, 1);

for array = 1:N_a
    snapshots_vs_arrays{array} = make_snapshots(tx_vs_arrays{array}, size_fft, overlap, window);
    fprintf('  Array %d: Created %d snapshots\n', array, size(snapshots_vs_arrays{array}, 3));
    
    r_vs_arrays{array} = create_vs_csm_6element(snapshots_vs_arrays{array}, bin_index, delta, ...
                                        bin_freqs, rho_0, N_v, num_bins, c_0, mic_positions_per_vs{array});
    fprintf('  Array %d: Created %dx%dx%d cross-spectral matrix\n', ...
            array, size(r_vs_arrays{array}));
end

%% PERFORM MVDR BEAMFORMING %%

fprintf('\n<strong>Running MVDR Beamforming (6-element)...</strong>\n');
fprintf('  Regularisation parameter: %.2e\n', loading);

response_db_arrays = cell(N_a, 1);

fprintf('\nCalculating individual array responses:\n');
for array = 1:N_a
    fprintf('  Processing Array %d...\n', array);
    [response_db_arrays{array}, mean_condition_number] = mvdr_vs_beamforming_6element(r_vs_arrays{array}, ...
        vs_centres_arrays{array}, mic_positions_per_vs{array}, candidate_points, bin_freqs, c_0, rho_0, ...
        loading, num_bins, N_v);
    fprintf('    Array %d response range: [%.2f, %.2f] dB\n', array, ...
        min(response_db_arrays{array}), max(response_db_arrays{array}));
end

% 2D scans for individual arrays
fprintf('\n2D Scans - Individual Arrays:\n');
for array = 1:N_a
    fprintf('  Plotting Array %d response...\n', array);
    [~, ~] = plot_2dscan_mvdr_6element(r_vs_arrays{array}, vs_centres_arrays{array}, ...
        mic_positions_per_vs{array}, candidate_points, bin_freqs, c_0, rho_0, y_scan, x_scan, ...
        X_grid, Y_grid, loading, num_bins, N_v, num_sources, source_positions, ...
        response_db_arrays{array}, array, grid_res, source_frequencies(1));

    if exist('batch_save_figures', 'var') && batch_save_figures
        saveas(gcf, [batch_figure_prefix '_2d_scan_6elem.png']);
    end
end

if N_a > 1
    fprintf('\nCombining array responses...\n');

    response_linear_combined = 1;
    for k = 1:N_a
        response_linear_combined = response_linear_combined .* 10.^(response_db_arrays{k}/10);
    end

    response_db_combined = 10*log10(response_linear_combined + eps);
    response_db_combined = response_db_combined - max(response_db_combined);

    fprintf('  Combined response range: [%.2f, %.2f] dB\n', ...
        min(response_db_combined), max(response_db_combined));

    fprintf('  Plotting superimposed response and localising sources...\n');
    [est_positions_combined, scan_2d_max_db] = ...
        superimposed_localisation(response_db_combined, X_grid, Y_grid, y_scan, x_scan, ...
        num_sources, source_positions, N_a, grid_res);

    if exist('batch_save_figures', 'var') && batch_save_figures
        saveas(fig_handle, [batch_figure_prefix '_2d_scan_combined_6elem.png']);
    end
end

%% PERFORMANCE METRICS %%

fprintf('\n<strong>PERFORMANCE METRICS</strong>\n');

if N_a == 1
    response_names = {'Array 1'};
    responses_to_analyse = {response_db_arrays{1}};
    angular_refs = {[vs_centres_arrays{1}(1:2, 1)]};
else
    response_names = [arrayfun(@(i) sprintf('Array %d', i), 1:N_a, 'UniformOutput', false), 'Superimposed'];
    responses_to_analyse = [response_db_arrays', {response_db_combined}];
    
    angular_refs = [cellfun(@(vs) vs(1:2, 1), vs_centres_arrays, 'UniformOutput', false)', ...
                   {mean(array_centres(1:2, :), 2)}];
end

num_cases = length(responses_to_analyse);
all_cases_results = cell(num_cases, 1);

min_separation = max([0.25, 3*grid_res]);
if num_sources > 1
    fprintf('Minimum source separation for multi-source detection: %.3f m\n', min_separation);
end

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
    response_work = grid_response;
    
    for src = 1:num_sources
        [peak_responses(src), max_idx] = max(response_work(:));
        [y_idx, x_idx] = ind2sub(size(grid_response), max_idx);
        est_positions(src, :) = [X_grid(y_idx, x_idx), Y_grid(y_idx, x_idx)];
        
        if num_sources > 1 && src < num_sources
            distances = sqrt((X_grid - est_positions(src,1)).^2 + (Y_grid - est_positions(src,2)).^2);
            response_work(distances < min_separation) = -Inf;
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
        fprintf('  Mean Squared Error: %.4f m²\n', MSE);
        fprintf('  MSE (normalised): %.3f%% of grid resolution\n', MSE_percent);
        fprintf('  Angular Error: %.3f deg\n', angular_error_deg);
        fprintf('  Peak Response: %.3f dB\n\n', peak_responses(src));
    end
    
    all_cases_results{case_idx} = results;
end

if N_a > 1
    fprintf('\n<strong>COMPARATIVE ANALYSIS:</strong>\n\n');
    
    print_comparison_table = @(metric_name, format_spec, get_metric) ...
        fprintf(['%s Comparison:\n%-10s ', repmat('%-15s ', 1, num_cases), '\n%s\n'], ...
                metric_name, 'Source', response_names{:}, repmat('-', 1, 10 + 15*num_cases)) && ...
        arrayfun(@(src) fprintf(['%-10d ', repmat(format_spec, 1, num_cases), '\n'], src, ...
                 arrayfun(@(c) get_metric(all_cases_results{c}, src), 1:num_cases)), 1:num_sources) && ...
        fprintf('\n');
    
    print_comparison_table('Radial Error (m)', '%-15.4f ', @(r, s) r.radial_error(s));
    print_comparison_table('Angular Error (deg)', '%-15.3f ', @(r, s) r.angular_error(s));
    
    fprintf('Performance Improvement:\n');
    fprintf('%-10s %-25s %-25s\n', 'Source', 'Radial Improvement', 'Angular Improvement');
    fprintf('%s\n', repmat('-', 1, 65));
    
    for src = 1:num_sources
        individual_radial = arrayfun(@(i) all_cases_results{i}.radial_error(src), 1:N_a);
        individual_angular = arrayfun(@(i) all_cases_results{i}.angular_error(src), 1:N_a);
        best_radial = min(individual_radial);
        best_angular = min(individual_angular);
        
        super_radial = all_cases_results{end}.radial_error(src);
        super_angular = all_cases_results{end}.angular_error(src);
        
        radial_improvement = best_radial - super_radial;
        radial_pct = 100 * radial_improvement / best_radial;
        angular_improvement = best_angular - super_angular;
        angular_pct = 100 * angular_improvement / best_angular;
        
        fprintf('%-10d %-25s %-25s\n', src, ...
            sprintf('%.4f m (%.1f%%)', radial_improvement, radial_pct), ...
            sprintf('%.3f deg (%.1f%%)', angular_improvement, angular_pct));
    end
    fprintf('\n');
end

all_results = all_cases_results{end};


%% BEAM PATTERN ANALYSIS %%
fprintf('\n<strong>GENERATING BEAM PATTERNS (6-element)</strong>\n');

distances_all = [];
for array = 1:N_a
    array_centre_2d = vs_centres_arrays{array}(1:2, 1);
    distances = sqrt(sum((source_positions(:,1:2) - array_centre_2d').^2, 2));
    distances_all = [distances_all; distances];
end
radius = median(distances_all);
fprintf('Using radius = %.3f m (median source distance across all arrays)\n', radius);

for array = 1:N_a
    array_centre_2d = vs_centres_arrays{array}(1:2, 1);
    [theta_deg, beam_pattern, fig_handle, beam_metrics] = compute_nearfield_beam_pattern_6element(...
        r_vs_arrays{array}, vs_centres_arrays{array}, mic_positions_per_vs{array}, array_centre_2d, ...
        source_positions(:,1:2), bin_freqs, c_0, rho_0, loading, num_bins, N_v, ...
        radius, sprintf('Array %d Near-Field Beam Pattern (6-elem, f = %.0f Hz)', array, source_frequencies(1)));

    if exist('batch_save_figures', 'var') && batch_save_figures
        saveas(fig_handle, [batch_figure_prefix '_beam_pattern_6elem.png']);
    end
end

if N_a > 1
    fprintf('\nComputing beam pattern for superimposed arrays (6-element)...\n');
    
    reference_centre = mean(array_centres(1:2, :), 2);
    
    r_vs_combined = zeros(size(r_vs_arrays{1}));
    for array = 1:N_a
        r_vs_combined = r_vs_combined + r_vs_arrays{array};
    end
    r_vs_combined = r_vs_combined / N_a;
    
    vs_centres_combined = [];
    mic_positions_combined = [];
    for array = 1:N_a
        vs_centres_combined = [vs_centres_combined, vs_centres_arrays{array}];
        for vs = 1:N_v
            mic_positions_combined = cat(3, mic_positions_combined, mic_positions_per_vs{array}(:,:,vs));
        end
    end
    
    N_v_combined = N_v * N_a;
    [theta_deg, beam_pattern, fig_handle, beam_metrics] = compute_nearfield_beam_pattern_6element(...
        r_vs_combined, vs_centres_combined, mic_positions_combined, reference_centre, ...
        source_positions(:,1:2), bin_freqs, c_0, rho_0, loading, num_bins, ...
        N_v_combined, radius, 'Superimposed Near-Field Beam Pattern (6-elem)');

    if exist('batch_save_figures', 'var') && batch_save_figures
        saveas(fig_handle, [batch_figure_prefix '_beam_pattern_superimposed_6elem.png']);
    end
end

%% RESTORE FIGURE VISIBILITY %%

if batch_mode
    set(0, 'DefaultFigureVisible', 'on');
end


%% FUNCTION DEFINITIONS %%

% FUNCTION: MAKE SNAPSHOTS OF GENERATED SIGNAL
function snapshots = make_snapshots(tx_vs, size_fft, overlap, window)
    step = round(overlap * size_fft);
    num_snap = floor((size(tx_vs,2) - size_fft) / step) + 1;
    start_idx = 1 + (0:(num_snap-1)) * step;
    
    idx_matrix = start_idx + (0:size_fft-1)';
    snapshots = tx_vs(:, idx_matrix);
    snapshots = reshape(snapshots, size(tx_vs,1), size_fft, num_snap);
    snapshots = snapshots .* reshape(window, [1, size_fft, 1]);
end


% FUNCTION: CREATE CROSS-SPECTRAL MATRIX FOR 6-ELEMENT VECTOR SENSOR
function r_vs = create_vs_csm_6element(snapshots_vs, bin_index, delta, bin_freqs, rho_0, N_v, num_bins, c_0, mic_positions_vs)
    signal_dim = 6 * N_v;
    
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
            
            [p0_vs, p1_vs, p2_vs, p3_vs, vx_vs, vy_vs] = vs_outputs_6element(fx_snap, delta, freq, rho_0, N_v, c_0);
            
            rho_c = rho_0 * c_0;
            vx_vs = rho_c * vx_vs;
            vy_vs = rho_c * vy_vs;
            
            vs_signal = zeros(signal_dim, 1);
            for n = 1:N_v
                idx = (n-1)*6 + (1:6);
                vs_signal(idx) = [p0_vs(n); p1_vs(n); p2_vs(n); p3_vs(n); vx_vs(n); vy_vs(n)];
            end
            
            R_freq = R_freq + (vs_signal * vs_signal');
        end
        
        r_vs(:,:,jf) = R_freq / num_snap;
    end
end


% FUNCTION: 6-ELEMENT VECTOR SENSOR OUTPUTS
function [p0_vs, p1_vs, p2_vs, p3_vs, vx_vs, vy_vs] = vs_outputs_6element(tx_freq, delta, freq, rho_0, N_v, c_0)
    omega = 2*pi*freq;
    
    p0_vs = zeros(N_v, 1);
    p1_vs = zeros(N_v, 1);
    p2_vs = zeros(N_v, 1);
    p3_vs = zeros(N_v, 1);
    vx_vs = zeros(N_v, 1);
    vy_vs = zeros(N_v, 1);
    
    for vs = 1:N_v
        idx = (vs-1)*4 + (1:4);
        
        p0_vs(vs) = tx_freq(idx(1));
        p1_vs(vs) = tx_freq(idx(2));
        p2_vs(vs) = tx_freq(idx(3));
        p3_vs(vs) = tx_freq(idx(4));
        
        mic_positions = [
            -delta/2, -delta/2;
             delta/2, -delta/2;
             delta/2,  delta/2;
            -delta/2,  delta/2
        ];
        
        M = [ones(4, 1), mic_positions];
        p_vector = [p0_vs(vs); p1_vs(vs); p2_vs(vs); p3_vs(vs)];
        
        coeffs = M \ p_vector;
        dpdx = coeffs(2);
        dpdy = coeffs(3);
        
        vx_vs(vs) = -dpdx / (1i * omega * rho_0);
        vy_vs(vs) = -dpdy / (1i * omega * rho_0);
    end
end


% FUNCTION: MVDR BEAMFORMING FOR 6-ELEMENT VECTOR SENSOR ARRAY
function [response_db, mean_condition_number] = mvdr_vs_beamforming_6element(r_vs, vs_centres, mic_positions_vs, candidate_points, ...
    bin_freqs, c_0, rho_0, loading, num_bins, N_v)
    
    num_points = size(candidate_points, 1);
    mvdr_responses = zeros(num_points, num_bins);
    condition_numbers = zeros(num_bins, 1);
    
    for jf = 1:num_bins 
        freq = bin_freqs(jf);
        rf = squeeze(r_vs(:,:,jf));
        
        [ur, er] = eig(rf);
        erv = real(diag(er));
        erv = max(erv, 1e-12);
        condition_numbers(jf) = max(erv) / min(erv(erv > 0));

        for n = 1:num_points
            source_pos = candidate_points(n,:).';
            
            v_vs = vs_steering_vector_6element(vs_centres, mic_positions_vs, source_pos, freq, c_0, rho_0, N_v);
            
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
    peak_response = max(responses_sum);
    if peak_response > 0
        responses_sum = max(responses_sum, peak_response * 1e-10);
    end
    response_db = 10*log10(responses_sum + eps);
    response_db = response_db - max(response_db);

    if nargout > 1
        mean_condition_number = mean(condition_numbers);
    else
        mean_condition_number = NaN;
    end
end


% FUNCTION: 6-ELEMENT STEERING VECTOR FOR ACOUSTIC VECTOR ARRAY
function v_vs = vs_steering_vector_6element(vs_centres, mic_positions_vs, source_pos, freq, c_0, rho_0, N_v)
    k = 2*pi*freq/c_0;
    omega = 2*pi*freq;
    
    v_vs = zeros(6*N_v, 1);
    
    for n = 1:N_v
        mic_pos = mic_positions_vs(:, :, n);
        
        p_steer = zeros(4, 1);
        for mic = 1:4
            r_vec = mic_pos(:, mic) - source_pos;
            r = norm(r_vec);
            p_steer(mic) = exp(-1i*k*r) / r;
        end
        
        r_vec_center = vs_centres(:,n) - source_pos;
        r_center = norm(r_vec_center);
        r_hat = r_vec_center / r_center;
        
        common_term = exp(-1i*k*r_center) / (1i*omega*rho_0*r_center^2);
        velocity_factor = (1 + 1i*k*r_center);
        
        vx_steer = common_term * velocity_factor * r_hat(1);
        vy_steer = common_term * velocity_factor * r_hat(2);
        
        rho_c = rho_0*c_0;
        vx_steer = rho_c * vx_steer;
        vy_steer = rho_c * vy_steer;
        
        idx = (n-1)*6 + (1:6);
        v_vs(idx) = [p_steer(1); p_steer(2); p_steer(3); p_steer(4); vx_steer; vy_steer];
    end
end


% FUNCTION: 2D SCAN PLOTTING FOR 6-ELEMENT MVDR
function [all_est_positions, scan_2d_max_db] = plot_2dscan_mvdr_6element(...
    r_vs, vs_centres, mic_positions_vs, candidate_points, bin_freqs, c_0, rho_0, y_scan, x_scan, ...
    X_grid, Y_grid, loading, num_bins, N_v, num_sources, source_positions, response_db, ...
    array_num, grid_res, source_freq)
    
    grid_response = reshape(response_db, length(y_scan), length(x_scan));
    
    figure('Name', sprintf('Array %d: 2D MVDR Beamforming (6-element)', array_num));
    imagesc(x_scan, y_scan, grid_response);
    axis xy;
    xlabel('X (m)', 'FontName', 'Times New Roman', 'FontSize', 14);
    ylabel('Y (m)', 'FontName', 'Times New Roman', 'FontSize', 14);
    colorbar;
    colormap('jet');
    hold on;
    
    [scan_2d_max_db, ~] = max(grid_response(:));
    
    min_separation = max([0.25, 3*grid_res]);
    
    all_est_positions = zeros(num_sources, 2);
    response_work = grid_response;
    
    for src = 1:num_sources
        [~, max_idx] = max(response_work(:));
        [y_idx, x_idx] = ind2sub(size(grid_response), max_idx);
        all_est_positions(src, :) = [X_grid(y_idx, x_idx), Y_grid(y_idx, x_idx)];
        
        if num_sources > 1 && src < num_sources
            distances = sqrt((X_grid - all_est_positions(src,1)).^2 + (Y_grid - all_est_positions(src,2)).^2);
            response_work(distances < min_separation) = -Inf;
        end
    end
    
    source_colours = lines(num_sources);
    
    for src = 1:num_sources
        true_x = source_positions(src, 1);
        true_y = source_positions(src, 2);
        est_x_src = all_est_positions(src, 1);
        est_y_src = all_est_positions(src, 2);
        
        plot(true_x, true_y, 'x', 'Color', source_colours(src,:), ...
            'MarkerSize', 14, 'LineWidth', 3, ...
            'DisplayName', sprintf('True Source %d', src));
        
        plot(est_x_src, est_y_src, 'o', 'Color', source_colours(src,:), ...
            'MarkerSize', 12, 'LineWidth', 2, 'MarkerFaceColor', 'none', ...
            'DisplayName', sprintf('Est Source %d', src));
    end
    
    legend('Location', 'northeast', 'Color', 'w', 'TextColor', 'k');
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
    
    title(sprintf('Array %d: MVDR Beamforming 6-elem (%d Sources, f = %.0f Hz)', ...
        array_num, num_sources, source_freq), ...
        'FontName', 'Times New Roman', 'FontSize', 14);
end


% FUNCTION: SUPERIMPOSED LOCALISATION
function [all_est_positions, max_db] = ...
    superimposed_localisation(response_db_combined, X_grid, Y_grid, y_scan, x_scan, ...
    num_sources, source_positions, N_a, grid_res)
    
    grid_response = reshape(response_db_combined, length(y_scan), length(x_scan));
    
    [max_db, ~] = max(grid_response(:));
    
    min_separation = max([0.25, 3*grid_res]);
    
    all_est_positions = zeros(num_sources, 2);
    response_work = grid_response;
    
    for src = 1:num_sources
        [~, max_idx] = max(response_work(:));
        [y_idx, x_idx] = ind2sub(size(grid_response), max_idx);
        all_est_positions(src, :) = [X_grid(y_idx, x_idx), Y_grid(y_idx, x_idx)];
        
        if num_sources > 1 && src < num_sources
            distances = sqrt((X_grid - all_est_positions(src,1)).^2 + (Y_grid - all_est_positions(src,2)).^2);
            response_work(distances < min_separation) = -Inf;
        end
    end
    
    figure('Name', sprintf('Superimposed MVDR (6-elem, %d Arrays)', N_a));
    imagesc(x_scan, y_scan, grid_response);
    axis xy;
    xlabel('X (m)', 'FontName', 'Times New Roman', 'FontSize', 14);
    ylabel('Y (m)', 'FontName', 'Times New Roman', 'FontSize', 14);
    colorbar;
    colormap('jet');
    hold on;
    
    source_colours = lines(num_sources);
    
    for src = 1:num_sources
        true_x = source_positions(src, 1);
        true_y = source_positions(src, 2);
        est_x_src = all_est_positions(src, 1);
        est_y_src = all_est_positions(src, 2);
        
        plot(true_x, true_y, 'x', 'Color', source_colours(src, :), ...
            'MarkerSize', 14, 'LineWidth', 3, ...
            'DisplayName', sprintf('True Source %d', src));
        
        plot(est_x_src, est_y_src, 'o', 'Color', source_colours(src, :), ...
            'MarkerSize', 12, 'LineWidth', 2, 'MarkerFaceColor', 'none', ...
            'DisplayName', sprintf('Est Source %d', src));
        
        fprintf('    Source %d: Estimated at (%.4f, %.4f) m\n', src, est_x_src, est_y_src);
    end
    
    legend('Location', 'northeast', 'Color', 'w', 'TextColor', 'k');
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
    title(sprintf('Superimposed Response (6-elem, %d Arrays) - Source Localisation', N_a), ...
        'FontName', 'Times New Roman', 'FontSize', 14);
end


% FUNCTION: COMPUTE NEAR-FIELD BEAM PATTERN (6-element version)
function [theta_deg, beam_pattern, fig_handle, beam_metrics] = compute_nearfield_beam_pattern_6element(...
    r_vs, vs_centres, mic_positions_vs, array_centre, source_positions, bin_freqs, c_0, rho_0, ...
    loading, num_bins, N_v, radius, plot_title)
    
    fprintf('\n<strong>Computing Near-Field Beam Pattern (6-element):</strong>\n');
    fprintf('  Range: %.3f m\n', radius);
    fprintf('  Array centre: (%.3f, %.3f) m\n', array_centre);
    
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
        
        lambda = loading * max(erv);
        
        for angle_idx = 1:num_angles
            x_pos = array_centre(1) + radius * sin(theta_rad(angle_idx));
            y_pos = array_centre(2) + radius * cos(theta_rad(angle_idx));
            source_pos = [x_pos; y_pos; 0];
            
            v_vs = vs_steering_vector_6element(vs_centres, mic_positions_vs, source_pos, freq, c_0, rho_0, N_v);
            
            rxv = (ur * diag(1./(erv + lambda)) * ur') * v_vs;
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
    peak_linear = max(beam_pattern_linear);
    if peak_linear > 0
        beam_pattern_linear = max(beam_pattern_linear, peak_linear * 1e-10);
    end
    beam_pattern = 10*log10(beam_pattern_linear + eps);
    beam_pattern = beam_pattern - max(beam_pattern);
    
    % Create polar plot
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
    
    % Calculate and store beam metrics
    beam_metrics = struct();
    
    [~, peak_idx] = max(beam_pattern);
    fprintf('\n<strong>Beam Pattern Statistics:</strong>\n');
    fprintf('  Peak response: 0.00 dB at %.1f degrees\n', theta_deg(peak_idx));
    
    % Calculate -3 dB beamwidth
    threshold_3db = max(beam_pattern) - 3;
    above_threshold = beam_pattern >= threshold_3db;
    above_threshold = above_threshold(:);
    diff_threshold = diff([0; above_threshold; 0]);
    starts = find(diff_threshold == 1);
    ends = find(diff_threshold == -1) - 1;
    
    beam_metrics.beamwidth_3dB = NaN;
    
    if ~isempty(starts)
        fprintf('  -3 dB Beamwidth(s):\n');
        for i = 1:length(starts)
            if ends(i) >= starts(i)
                beamwidth = theta_deg(ends(i)) - theta_deg(starts(i));
                if beamwidth < 0
                    beamwidth = beamwidth + 360;
                end
                if i == 1
                    beam_metrics.beamwidth_3dB = beamwidth;
                end
                fprintf('    Region %d: %.1f degrees (%.1f° to %.1f°)\n', ...
                    i, beamwidth, theta_deg(starts(i)), theta_deg(ends(i)));
            end
        end
    end
    
    % Calculate sidelobe level
    main_lobe_mask = beam_pattern >= threshold_3db;
    main_lobe_mask = main_lobe_mask(:);
    beam_metrics.sidelobe_level = NaN;
    if any(~main_lobe_mask)
        beam_pattern_col = beam_pattern(:);
        sidelobe_level = max(beam_pattern_col(~main_lobe_mask));
        beam_metrics.sidelobe_level = sidelobe_level;
        fprintf('  Maximum sidelobe level: %.2f dB\n', sidelobe_level);
    end
    
    % Calculate directivity index
    beam_linear = 10.^(beam_pattern(:)/10);
    theta_rad_col = theta_rad(:);
    integral_val = trapz(theta_rad_col, beam_linear);
    DI = 10*log10(2*pi / integral_val);
    beam_metrics.directivity_index = DI;
    fprintf('  Directivity Index: %.2f dB\n', DI);
    
    fprintf('\n');
end