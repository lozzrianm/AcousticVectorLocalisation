clc; clear all; close all;

% L Marshall 27/11/2025
% Variable Number of Sources Localisation via an Acoustic Vector Array

% MVDR method of beamforming
% Using broadband analysis with spectrum blocks
% with eigenvalue decomposition, adaptive regularisation, and weighted MVDR

%% DEFINE INPUT VARIABLES %%

csv_filename = 'generatedsignal_avs.csv';

% Physical Environment Parameters
c_0 = 340; %Speed of sound (m/s)
rho_0 = 1.02; %Air density at STP

% Visualisation Parameters
% Source positions [x, y, z] (m)
source_positions = [
    0, -0.25, 0; %Source 1
    %4, -2, 0.0; %Source 2
];

% Source frequencies (Hz)
source_frequencies = [
    1000; %Source 1 freq
    %2000 %Source 2 freq
];

% Source amplitudes 
source_amplitudes = [
    0.01; %Source 1 a
    %0.01; %Source 2 a
];

num_sources = size(source_positions, 1);

% Vector Sensor Array(s) Parameters
N_a = 1; %Number of independent arrays
N_v = 1; %No. of vector sensors per array
d_y = 0.1; %Vector sensor y axis spacing (m)
delta = 0.0425; %MEMS colocation spacing (m)

% Array centre geom 1
x_a1 = 2.0 + delta/2; % Geometric centre x-coordinate
z_a1 = 0; % Array centre z-coordinate
y_a1 = -d_y/2 + delta/2; % Geometric centre y-coordinate
A1_coord = [x_a1; y_a1; z_a1];

% Array centre geom 2
lambda = c_0/source_frequencies(1,:);
a_spacing = sqrt(((10*lambda)^2)/2); %x and y spacing for 10*lambda distance
x_a2 = x_a1 - a_spacing; % Geometric centre x-coordinate
z_a2 = 0; % Array centre z-coordinate
y_a2 = y_a1 - a_spacing; % Geometric centre y-coordinate
A2_coord = [x_a2; y_a2; z_a2];

% Sampling Characteristics
freq_limit = 3000; %Max freq for analysis (Hz)
overlap = 0.5; %Percentage overlap for sample window (%)
loading = 1e-4; %Regularisation parameter

% Grid Search Parameters
x_scan_points = 200; 
y_scan_points = 200;
x_margin = 3; %(m)
y_margin = 2; %(m)

%% LOAD GENERATED SIGNAL FROM CSV %%
% Data generated in 'VA_Output_WIP.m'

fprintf('Loading signal data from: %s\n', csv_filename);

data_table = readtable(csv_filename);
data = table2array(data_table);
time = data(:,1); %Time vector
[numRows, numCols] = size(data);

N = numRows; %Number of samples
N_ma = N_v * 4; %Microphones per array
N_m = N_a * N_ma; %Total microphones across both arrays

fprintf('Loaded %d samples from %d microphones (%d arrays, %d mics each)\n', ...
    N, N_m, N_a, N_ma);

%% COMBINE SIGNAL COMPONENTS FOR PROCESSING %%
% Combine monopole and sinusoidal signals for each array separately

% Preallocate cell array
tx_vs_arrays = cell(N_a, 1);

% Process signals
for array = 1:N_a
    tx_vs = zeros(N_ma, N);
    
    for idx_local = 1:N_ma
        % Calculate column indices in CSV
        monopole_col = 1 + (array-1)*N_ma + idx_local;
        sinusoidal_col = 1 + N_m + (array-1)*N_ma + idx_local;
        % Combine monopole and sinusoidal components
        tx_vs(idx_local, :) = data(:, monopole_col) + data(:, sinusoidal_col);
    end
    tx_vs_arrays{array} = tx_vs;
end

%% ARRAY GEOMETRY SET UP %%
% Combine array centres 
array_centres = [A1_coord, A2_coord];

% Adjust for if multiple VS per array
if N_v > 1
    array_centres(2, :) = array_centres(2, :) - (N_v-1)*d_y/2;
end

vs_centres_arrays = cell(N_a, 1);
mic_positions = zeros(3, N_m);
% To define 2x2 square config of mics in x,y,z
mic_offsets = delta/2 * [-1, 1, 1, -1; -1, -1, 1, 1; 0, 0, 0, 0];

for array = 1:N_a
    vs_centres = zeros(3, N_v);
    
    for vs = 1:N_v
        % VS centre
        vs_centres(:, vs) = array_centres(:, array) + [0; (vs-1)*d_y; 0];
        
        % Microphone positions
        idx = (array-1)*N_ma + (vs-1)*4 + (1:4);
        mic_positions(:, idx) = vs_centres(:, vs) + mic_offsets;
    end
    
    vs_centres_arrays{array} = vs_centres;
end

if N_a > 1
    fprintf('Array separation: %.3f m\n', norm(diff(array_centres, 1, 2)));
end

%% CONVERT TO FREQUENCY DOMAIN %%

fprintf('Converting to frequency domain...\n');

d_t = time(2) - time(1); %Time step (s)
F_s = 1/d_t; %Sampling frequency (Hz)
f = F_s * (0:floor(N/2)) / N; %Frequency vector

% Limit frequency range
idx_limit = f <= freq_limit;

% Store time-domain signals in cell array for looping
tx_vs_arrays_time = tx_vs_arrays;

% Preallocate FFT storage
tx_freq_arrays = cell(N_a, 1);

% Create figure for frequency spectra
figure('Name', 'Frequency Domain - All Arrays');

% Process each array
for array = 1:N_a
    % Preallocate FFT matrix for this array
    tx_freq_arrays{array} = zeros(N_ma, sum(idx_limit));
    
    % Create subplot for this array
    subplot(N_a, 1, array);
    hold on;
    
    % Process each sensor in this array
    for sensor = 1:N_ma
        % FFT of each sensor signal
        Y = fft(tx_vs_arrays_time{array}(sensor, :));
        Y = Y(1:floor(N/2)+1); %Keep positive freqs (complex)
        
        % Store only the useful freqs (<= freq_limit)
        tx_freq_arrays{array}(sensor, :) = Y(idx_limit);
        
        % Plot magnitude spectrum
        Y_mag = abs(Y)/N;
        Y_mag(2:end-1) = 2*Y_mag(2:end-1);
        Y_dB = 20*log10(Y_mag + eps);
        plot(f, Y_dB, 'DisplayName', sprintf('Sensor %d', sensor));
    end
    
    hold off;
    xlabel('Frequency (Hz)', 'FontName', 'Times New Roman', 'FontSize', 14);
    ylabel('Amplitude (dB)', 'FontName', 'Times New Roman', 'FontSize', 14);
    title(sprintf('Array %d Frequency Spectrum', array), ...
          'FontName', 'Times New Roman', 'FontSize', 14);
    xlim([0 freq_limit]);
    legend('show', 'Location', 'best');
    grid on;
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
end

%% DEFINING THE GRID SEARCH AREA %%

fprintf('Defining search area...\n');

% Calculate centroid of arrays for grid centring
x_arrays_mean = mean(array_centres(1, :));
y_arrays_mean = mean(array_centres(2, :));

if N_a == 1
    x_scan = linspace(x_a1 - x_margin, x_a1 + x_margin, x_scan_points);
    y_scan = linspace(y_a1 - y_margin, y_a1 + y_margin, y_scan_points);
    z_scan = 0; %2D plane, change if needed
else
% Use array mean as grid centre
    x_scan = linspace(x_arrays_mean - x_margin, x_arrays_mean + x_margin, x_scan_points);
    y_scan = linspace(y_arrays_mean - y_margin, y_arrays_mean + y_margin, y_scan_points);
    z_scan = 0; %2D plane, change if needed
end

[X_grid, Y_grid] = meshgrid(x_scan, y_scan);
Z_grid = zeros(size(X_grid));

% N_m x 3 matrix
candidate_points = [X_grid(:), Y_grid(:), Z_grid(:)];

% Grid resolution in x and y directions
x_res = x_scan(2) - x_scan(1);
y_res = y_scan(2) - y_scan(1);
grid_res = mean([x_res, y_res]);

fprintf('Search grid: %d x %d points, resolution: %.4f m\n', ...
    length(x_scan), length(y_scan), grid_res);

%% CREATING FREQUENCY BINS FOR BROADBAND ANALYSIS %%

fprintf('Creating frequency bins...\n');

% Maximum bin width for required spatial resolution
max_binwidth = 1/(8*d_y*(N_ma-1)/c_0); %(Hz)
size_fft = floor(F_s/max_binwidth);

% Frequency spacing
delta_f = F_s / N; %Frequency spacing per bin
fft_vec = F_s * (0:(size_fft-1)) / size_fft; %Freq vector

% Target frequencies
% Step in terms of FFT bins
min_freq = max_binwidth;
target_freqs = min_freq:max_binwidth:freq_limit;

% Map each target to nearest snapshot FFT bin
bin_index = zeros(size(target_freqs));
for k = 1:length(target_freqs)
    [~, bin_index(k)] = min(abs(fft_vec - target_freqs(k)));
end

% Generate bin indices within accepted range
bin_index = unique(bin_index); %Remove duplicates
bin_index = bin_index(bin_index <= size_fft);
bin_freqs = fft_vec(bin_index); %Actual bin freqs
num_bins = length(bin_index); %No. bins

fprintf('Using %d frequency bins (%.1f - %.1f Hz)\n', ...
    num_bins, min(bin_freqs), max(bin_freqs));

num_snaps = 8*N_m; %No. snapshots

% Calc no. samples accounting for overlap
num_samps = (num_snaps+1)*size_fft*overlap;

% Taper window
window = hanning(size_fft)';

%% CREATE AVS SNAPSHOTS AND MAKE CSM %%

fprintf('Creating snapshots and cross-spectral matrix...\n');

snapshots_vs_arrays = cell(N_a, 1);
r_vs_arrays = cell(N_a, 1);

for array = 1:N_a
    % Create snapshots
    snapshots_vs_arrays{array} = make_snapshots(tx_vs_arrays{array}, size_fft, overlap, window);
    fprintf('  Array %d: Created %d snapshots\n', array, size(snapshots_vs_arrays{array}, 3));
    
    % Create cross-spectral matrix
    r_vs_arrays{array} = create_vs_csm(snapshots_vs_arrays{array}, bin_index, delta, ...
                                        bin_freqs, rho_0, N_v, num_bins, c_0);
    fprintf('  Array %d: Created %dx%dx%d cross-spectral matrix\n', ...
            array, size(r_vs_arrays{array}));
end

%% PERFORM MVDR BEAMFORMING %%

fprintf('\n<strong>Running MVDR Beamforming...</strong>\n');
fprintf('  Regularisation parameter: %.2e\n', loading);

% Calculate beamforming response for each array separately
response_db_arrays = cell(N_a, 1);

fprintf('\nCalculating individual array responses:\n');
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
    [~, ~] = plot_2dscan_mvdr(r_vs_arrays{array}, vs_centres_arrays{array}, ...
        candidate_points, bin_freqs, c_0, rho_0, y_scan, x_scan, ...
        X_grid, Y_grid, loading, num_bins, N_v, num_sources, source_positions, ...
        response_db_arrays{array}, array);
end

if N_a > 1
    % Superimpose responses via multiplication (intersection method)
    fprintf('\nCombining array responses...\n');

    % Convert from dB to linear and multiply responses
    response_linear_combined = 1; %Identity for multiplication
    for k = 1:N_a
        response_linear_combined = response_linear_combined .* 10.^(response_db_arrays{k}/10);
    end

    % Convert back to dB and normalise
    response_db_combined = 10*log10(response_linear_combined + eps);
    response_db_combined = response_db_combined - max(response_db_combined);

    fprintf('  Combined response range: [%.2f, %.2f] dB\n', ...
        min(response_db_combined), max(response_db_combined));

    % Plot combined response and find source positions
    fprintf('  Plotting superimposed response and localising sources...\n');
    [est_positions_combined, scan_2d_max_db] = ...
        superimposed_localisation(response_db_combined, X_grid, Y_grid, y_scan, x_scan, ...
        num_sources, source_positions, N_a);
end

%% PERFORMANCE METRICS %%

fprintf('\n<strong>PERFORMANCE METRICS</strong>\n');

% Determine number of cases based on array configuration
if N_a == 1
    % Single array: only one response to analyse
    response_names = {'Array 1'};
    responses_to_analyse = {response_db_arrays{1}};
    angular_refs = {[vs_centres_arrays{1}(1:2, 1)]};
else
    % Multiple arrays: analyse each array plus superimposed
    response_names = [arrayfun(@(i) sprintf('Array %d', i), 1:N_a, 'UniformOutput', false), 'Superimposed'];
    responses_to_analyse = [response_db_arrays', {response_db_combined}];
    
    % Angular reference points: individual array centres + midpoint
    angular_refs = [cellfun(@(vs) vs(1:2, 1), vs_centres_arrays, 'UniformOutput', false)', ...
                   {mean(array_centres(1:2, :), 2)}];
end

num_cases = length(responses_to_analyse);
all_cases_results = cell(num_cases, 1);

% Analyse each case
for case_idx = 1:num_cases
    fprintf('\n<strong>%s:</strong>\n', response_names{case_idx});
    
    % Prepare grid and reference
    grid_response = reshape(responses_to_analyse{case_idx}, length(y_scan), length(x_scan));
    ref_pos = angular_refs{case_idx};
    
    % Initialise results structure
    results = struct('est_positions', zeros(num_sources, 2), ...
                     'true_positions', zeros(num_sources, 2), ...
                     'errors', zeros(num_sources, 2), ...
                     'MSE', zeros(num_sources, 1), ...
                     'MSE_percent', zeros(num_sources, 1), ...
                     'radial_error', zeros(num_sources, 1), ...
                     'angular_error', zeros(num_sources, 1), ...
                     'peak_response', zeros(num_sources, 1));
    
    % Process each source
    for src = 1:num_sources
        fprintf('  <strong>SOURCE %d: (%.3f, %.3f) m @ %.0f Hz</strong>\n', ...
            src, source_positions(src, 1:2), source_frequencies(src));
        
        % Define search region
        true_pos = source_positions(src, 1:2);
        search_radius = [0.5, 0.5]; % [x, y] search radius (m)
        
        % Create search mask
        region_mask = (X_grid >= true_pos(1) - search_radius(1)) & ...
                     (X_grid <= true_pos(1) + search_radius(1)) & ...
                     (Y_grid >= true_pos(2) - search_radius(2)) & ...
                     (Y_grid <= true_pos(2) + search_radius(2));
        
        % Find peak in search region
        region_response = grid_response;
        region_response(~region_mask) = -Inf;
        [max_val, max_idx] = max(region_response(:));
        [y_idx, x_idx] = ind2sub(size(grid_response), max_idx);
        est_pos = [X_grid(y_idx, x_idx), Y_grid(y_idx, x_idx)];
        
        % Calculate errors
        errors = true_pos - est_pos;
        radial_error = norm(errors);
        MSE = sum(errors.^2);
        MSE_percent = 100 * MSE / (grid_res^2);
        
        % Angular error
        true_angle = atan2(true_pos(2) - ref_pos(2), true_pos(1) - ref_pos(1));
        est_angle = atan2(est_pos(2) - ref_pos(2), est_pos(1) - ref_pos(1));
        angular_error_deg = rad2deg(abs(true_angle - est_angle));
        if angular_error_deg > 180
            angular_error_deg = 360 - angular_error_deg;
        end
        
        % Store results
        results.est_positions(src, :) = est_pos;
        results.true_positions(src, :) = true_pos;
        results.errors(src, :) = errors;
        results.MSE(src) = MSE;
        results.MSE_percent(src) = MSE_percent;
        results.radial_error(src) = radial_error;
        results.angular_error(src) = angular_error_deg;
        results.peak_response(src) = max_val;
        
        % Display results
        fprintf('  Estimated position: (%.4f, %.4f) m\n', est_pos);
        fprintf('  Position error: dx = %.4f m, dy = %.4f m\n', errors);
        fprintf('  Radial error: %.4f m\n', radial_error);
        fprintf('  Mean Squared Error: %.4f m²\n', MSE);
        fprintf('  MSE (normalised): %.3f%% of grid resolution\n', MSE_percent);
        fprintf('  Angular Error: %.3f deg\n', angular_error_deg);
        fprintf('  Peak Response: %.3f dB\n\n', max_val);
    end
    
    all_cases_results{case_idx} = results;
end

% Comparative analysis (only for multiple arrays)
if N_a > 1
    fprintf('\n<strong>COMPARATIVE ANALYSIS:</strong>\n\n');
    
    % Helper function to print comparison table
    print_comparison_table = @(metric_name, format_spec, get_metric) ...
        fprintf(['%s Comparison:\n%-10s ', repmat('%-15s ', 1, num_cases), '\n%s\n'], ...
                metric_name, 'Source', response_names{:}, repmat('-', 1, 10 + 15*num_cases)) && ...
        arrayfun(@(src) fprintf(['%-10d ', repmat(format_spec, 1, num_cases), '\n'], src, ...
                 arrayfun(@(c) get_metric(all_cases_results{c}, src), 1:num_cases)), 1:num_sources) && ...
        fprintf('\n');
    
    % Print comparison tables
    print_comparison_table('Radial Error (m)', '%-15.4f ', @(r, s) r.radial_error(s));
    print_comparison_table('Angular Error (deg)', '%-15.3f ', @(r, s) r.angular_error(s));
    
    % Performance improvement analysis
    fprintf('Performance Improvement:\n');
    fprintf('%-10s %-25s %-25s\n', 'Source', 'Radial Improvement', 'Angular Improvement');
    fprintf('%s\n', repmat('-', 1, 65));
    
    for src = 1:num_sources
        % Best individual array performance
        individual_radial = arrayfun(@(i) all_cases_results{i}.radial_error(src), 1:N_a);
        individual_angular = arrayfun(@(i) all_cases_results{i}.angular_error(src), 1:N_a);
        best_radial = min(individual_radial);
        best_angular = min(individual_angular);
        
        % Superimposed performance
        super_radial = all_cases_results{end}.radial_error(src);
        super_angular = all_cases_results{end}.angular_error(src);
        
        % Calculate improvements
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

% Store results for backward compatibility
all_results = all_cases_results{end};

%% BEAM PATTERN ANALYSIS %%
fprintf('\n<strong>GENERATING BEAM PATTERNS</strong>\n');

% For individual arrays
for array = 1:N_a
    array_centre = vs_centres_arrays{array}(1:2, 1);
    plot_beam_pattern(response_db_arrays{array}, X_grid, Y_grid, ...
        array_centre, sprintf('Array %d Beam Pattern', array), source_positions(:,1:2));
end

if N_a > 1
    % For superimposed response
    reference_centre = mean(array_centres(1:2, :), 2);
    plot_beam_pattern(response_db_combined, X_grid, Y_grid, ...
        reference_centre, 'Superimposed Beam Pattern', source_positions(:,1:2));
end


%% FUNCTION DEFINITIONS %%

% FUNCTION: MAKE SNAPSHOTS OF GENERATED SIGNAL
% Make time domain snapshots of input signal
function snapshots = make_snapshots(tx_vs, size_fft, overlap, window)
    % No. samples to move per snapshot
    step = round(overlap * size_fft);
    % No. of full snapshots that fit in signal with given step size
    num_snap = floor((size(tx_vs,2) - size_fft) / step) + 1; %No. snapshots
    % Start indices for each snapshot
    start_idx = 1 + (0:(num_snap-1)) * step;
    
    % Create index matrix (size_fft x num_snap)
    idx_matrix = start_idx + (0:size_fft-1)'; %Each column = indices for a snapshot
    % Extract snapshots (force reshape into 3D)
    snapshots = tx_vs(:, idx_matrix);
    snapshots = reshape(snapshots, size(tx_vs,1), size_fft, num_snap);
    % Apply window along snapshots dimension
    snapshots = snapshots .* reshape(window, [1, size_fft, 1]);
end


% FUNCTION: CREATE CROSS-SPECTRAL MATRIX FOR VECTOR SENSOR ARRAY
% Convert 4-mic snapshots to pressure and velocity, then create CSM
function r_vs = create_vs_csm(snapshots_vs, bin_index, delta, bin_freqs, rho_0, N_v, num_bins, c_0)
    % Create CSM for AVS array
    % snapshots: N_m x size_fft x num_snap 
    signal_dim = 3 * N_v;
    % FFT of snapshots
    fx = fft(snapshots_vs, [], 2);
    fx = fx(:, bin_index, :);
    % Initialise CSM
    r_vs = zeros(signal_dim, signal_dim, num_bins);
    num_snap = size(snapshots_vs, 3);
    
    for jf = 1:num_bins
        freq = bin_freqs(jf);
        % Average over snapshots
        fx_bin = squeeze(fx(:, jf, :)); %M_mics x N_snap
        % Initialise CSM for this frequency bin
        R_freq = zeros(signal_dim, signal_dim);
        % Loop over snapshots
        for snap = 1:num_snap
            fx_snap = fx_bin(:, snap); % Pressure at each mic for this snapshot
            % Compute VS outputs for THIS snapshot
            [p_vs, vx_vs, vy_vs] = vs_outputs(fx_snap, delta, freq, rho_0, N_v, c_0);
            
            % Apply impedance scaling (plane wave approximation)
            rho_c = rho_0 * c_0;
            vx_vs = rho_c * vx_vs;
            vy_vs = rho_c * vy_vs;
            % Form signal vector for this snapshot
            vs_signal = zeros(signal_dim, 1);
            for n = 1:N_v
                idx = (n-1)*3 + (1:3);
                vs_signal(idx) = [p_vs(n); vx_vs(n); vy_vs(n)];
            end
            % Accumulate outer product
            R_freq = R_freq + (vs_signal * vs_signal');
        end
        % Average over snapshots
        r_vs(:,:,jf) = R_freq / num_snap;
    end
end


% FUNCTION: MVDR BEAMFORMING FOR VECTOR SENSOR ARRAY
function response_db = mvdr_vs_beamforming(r_vs, vs_centres, candidate_points, ...
    bin_freqs, c_0, rho_0, loading, num_bins, N_v)
    
    num_points = size(candidate_points, 1);
    mvdr_responses = zeros(num_points, num_bins);
    
    for jf = 1:num_bins 
        freq = bin_freqs(jf);
        rf = squeeze(r_vs(:,:,jf)); %Covariance matrix for this freq bin
        % Eigenvalue decomposition (improve stability)
        [ur, er] = eig(rf);
        erv = real(diag(er));
        % Safety check for negative eigenvalues 
        erv = max(erv, 1e-12);
        
        for n = 1:num_points
            source_pos = candidate_points(n,:).';
            % Create VS steering vector
            v_vs = vs_steering_vector(vs_centres, source_pos, freq, c_0, rho_0, N_v);
            
            % Adaptive regularisation
            lambda = loading * max(erv); %lambda = diagonal loading factor
            % Weighted MVDR using eigendecomposition
            rxv = (ur * diag(1./(erv + lambda)) * ur') * v_vs;
            denominator = v_vs' * rxv; %Normalise to get weight vector
            
            if abs(denominator) > 1e-12
                vmvdr = rxv / denominator;
                % Calc weighted output power
                mvdr_responses(n, jf) = abs(vmvdr' * rf * vmvdr);
            end
        end
    end
    % Sum across frequencies
    responses_sum = sum(mvdr_responses, 2);
    response_db = 10*log10(responses_sum + eps);
    response_db = response_db - max(response_db);
end


% FUNCTION: 2D SCAN PLOTTING FOR MVDR
function [all_est_positions, scan_2d_max_db] = plot_2dscan_mvdr(...
    r_vs, vs_centres, candidate_points, bin_freqs, c_0, rho_0, y_scan, x_scan, ...
    X_grid, Y_grid, loading, num_bins, N_v, num_sources, source_positions, response_db, array_num)
    
    % If response_db not provided, calculate it
    if nargin < 16
        response_db = mvdr_vs_beamforming(r_vs, vs_centres, candidate_points, ...
            bin_freqs, c_0, rho_0, loading, num_bins, N_v);
        array_num = 1; %Default to array 1
    elseif nargin < 17
        array_num = 1; %Default if not specified
    end
    
    % Reshape responses to a grid for plotting
    grid_response = reshape(response_db, length(y_scan), length(x_scan));
    
    figure('Name', sprintf('Array %d: 2D MVDR Beamforming', array_num));
    imagesc(x_scan, y_scan, grid_response);
    axis xy;
    xlabel('X (m)', 'FontName', 'Times New Roman', 'FontSize', 14);
    ylabel('Y (m)', 'FontName', 'Times New Roman', 'FontSize', 14);
    colorbar;
    colormap('jet');
    hold on;
    
    % Find global maximum 
    [scan_2d_max_db, max_idx] = max(grid_response(:));
    [y_idx, x_idx] = ind2sub(size(grid_response), max_idx);
    
    % Find estimated position for each source
    all_est_positions = zeros(num_sources, 2);
    search_radius_x = 0.5;
    search_radius_y = 0.5;
    
    source_colours = lines(num_sources);
    
    for src = 1:num_sources
        true_x = source_positions(src, 1);
        true_y = source_positions(src, 2);
        
        % Define search region around true source
        x_mask = (X_grid >= true_x - search_radius_x) & (X_grid <= true_x + search_radius_x);
        y_mask = (Y_grid >= true_y - search_radius_y) & (Y_grid <= true_y + search_radius_y);
        region_mask = x_mask & y_mask;
        
        % Find peak in region
        region_response = grid_response;
        region_response(~region_mask) = -Inf;
        [~, max_idx_src] = max(region_response(:));
        [y_idx_src, x_idx_src] = ind2sub(size(grid_response), max_idx_src);
        est_x_src = X_grid(y_idx_src, x_idx_src);
        est_y_src = Y_grid(y_idx_src, x_idx_src);
        
        all_est_positions(src,:) = [est_x_src, est_y_src];
        
        % Plot true source position
        plot(true_x, true_y, 'x', 'Color', source_colours(src,:), ...
            'MarkerSize', 14, 'LineWidth', 3, ...
            'DisplayName', sprintf('True Source %d', src));
        
        % Plot estimated position
        plot(est_x_src, est_y_src, 'o', 'Color', source_colours(src,:), ...
            'MarkerSize', 12, 'LineWidth', 2, 'MarkerFaceColor', 'none', ...
            'DisplayName', sprintf('Est Source %d', src));
    end
    
    legend('Location', 'northeast', 'Color', 'w', 'TextColor', 'k');
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
    title(sprintf('Array %d: MVDR Beamforming (%d Sources)', array_num, num_sources), ...
        'FontName', 'Times New Roman', 'FontSize', 14);
end


% FUNCTION: SUPERIMPOSED LOCALISATION
% Plot combined array response and localise sources
function [all_est_positions, max_db] = ...
    superimposed_localisation(response_db_combined, X_grid, Y_grid, y_scan, x_scan, ...
    num_sources, source_positions, N_a)
    
    % Reshape response to grid
    grid_response = reshape(response_db_combined, length(y_scan), length(x_scan));
    
    % Find global maximum
    [max_db, max_idx] = max(grid_response(:));
    [y_idx, x_idx] = ind2sub(size(grid_response), max_idx);
    
    % Find estimated position for each source
    all_est_positions = zeros(num_sources, 2);
    search_radius_x = 0.5; %(m)
    search_radius_y = 0.5; %(m)
    
    source_colours = lines(num_sources);
    
    % Create combined response plot
    figure('Name', sprintf('Superimposed MVDR (%d Arrays)', N_a));
    imagesc(x_scan, y_scan, grid_response);
    axis xy;
    xlabel('X (m)', 'FontName', 'Times New Roman', 'FontSize', 14);
    ylabel('Y (m)', 'FontName', 'Times New Roman', 'FontSize', 14);
    colorbar;
    colormap('jet');
    hold on;
    
    for src = 1:num_sources
        true_x = source_positions(src, 1);
        true_y = source_positions(src, 2);
        
        % Define search region around this source
        x_mask = (X_grid >= true_x - search_radius_x) & (X_grid <= true_x + search_radius_x);
        y_mask = (Y_grid >= true_y - search_radius_y) & (Y_grid <= true_y + search_radius_y);
        region_mask = x_mask & y_mask;
        
        % Find peak in region
        region_response = grid_response;
        region_response(~region_mask) = -Inf;
        [~, max_idx_src] = max(region_response(:));
        [y_idx_src, x_idx_src] = ind2sub(size(grid_response), max_idx_src);
        est_x_src = X_grid(y_idx_src, x_idx_src);
        est_y_src = Y_grid(y_idx_src, x_idx_src);
        
        all_est_positions(src, :) = [est_x_src, est_y_src];
        
        % Plot true source position
        plot(true_x, true_y, 'x', 'Color', source_colours(src, :), ...
            'MarkerSize', 14, 'LineWidth', 3, ...
            'DisplayName', sprintf('True Source %d', src));
        
        % Plot estimated position
        plot(est_x_src, est_y_src, 'o', 'Color', source_colours(src, :), ...
            'MarkerSize', 12, 'LineWidth', 2, 'MarkerFaceColor', 'none', ...
            'DisplayName', sprintf('Est Source %d', src));
        
        fprintf('    Source %d: Estimated at (%.4f, %.4f) m\n', src, est_x_src, est_y_src);
    end
    
    legend('Location', 'northeast', 'Color', 'w', 'TextColor', 'k');
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
    title(sprintf('Superimposed Response (%d Arrays) - Source Localisation', N_a), ...
        'FontName', 'Times New Roman', 'FontSize', 14);
end


% FUNCTION: VECTOR SENSOR OUTPUTS WITH LEAST SQUARES GRADIENT
% Compute vector sensor pressure - particle velocity conversions
function [p_vs, vx_vs, vy_vs] = vs_outputs(tx_freq, delta, freq, rho_0, N_v, c_0)
    omega = 2*pi*freq;
    % Initialise outputs
    p_vs = zeros(N_v, 1);
    vx_vs = zeros(N_v, 1);
    vy_vs = zeros(N_v, 1);
    
    for vs = 1:N_v
        idx = (vs-1)*4 + (1:4);
        % Get pressure at 4 microphone positions (square corners)
        p0 = tx_freq(idx(1)); %bottom left
        p1 = tx_freq(idx(2)); %bottom right
        p2 = tx_freq(idx(3)); %top right
        p3 = tx_freq(idx(4)); %top left
        % Microphone positions relative to centre
        mic_positions = [
            -delta/2, -delta/2;
             delta/2, -delta/2;
             delta/2,  delta/2;
            -delta/2,  delta/2
        ];
        
        % Construct design matrix M for least-squares fit
        % Each row: [1, x, y] for linear model p = a0 + a1*x + a2*y
        M = [ones(4, 1), mic_positions];
        % Pressure vector (complex)
        p_vector = [p0; p1; p2; p3];
        % Least-squares solution: M * [p_centre; dp/dx; dp/dy] = p_vector
        % Solution: coeffs = M \ p_vector
        coeffs = M \ p_vector;
        % Extract results
        p_vs(vs) = coeffs(1); %Pressure at centre
        dpdx = coeffs(2); %Pressure gradient in x
        dpdy = coeffs(3); %Pressure gradient in y
        
        % Convert pressure gradients to particle velocities with Euler eq
        vx_vs(vs) = -dpdx / (1i * omega * rho_0);
        vy_vs(vs) = -dpdy / (1i * omega * rho_0);
    end
end


% FUNCTION: STEERING VECTOR FOR ACOUSTIC VECTOR ARRAY
function v_vs = vs_steering_vector(vs_centres, source_pos, freq, c_0, rho_0, N_v)
    % Output: [p1; vx1; vy1; p2; vx2; vy2; ...]
    k = 2*pi*freq/c_0;
    omega = 2*pi*freq;
    
    v_vs = zeros(3*N_v, 1); %3 components per vector sensor
    
    for n = 1:N_v
        % Vector from source to VS centre
        r_vec = vs_centres(:,n) - source_pos;
        r = norm(r_vec);
        r_hat = r_vec / r;
        
        % Pressure steering component
        p_steer = exp(-1i*k*r) / r;
        
        % Particle velocity steering components
        common_term = exp(-1i*k*r) / (1i*omega*rho_0*r^2);
        velocity_factor = (1 + 1i*k*r);
        
        vx_steer = common_term * velocity_factor * r_hat(1);
        vy_steer = common_term * velocity_factor * r_hat(2);
        
        % Apply scaling factor to velocity terms 
        rho_c = rho_0*c_0;
        vx_steer = rho_c * vx_steer;
        vy_steer = rho_c * vy_steer;
        
        % Store in steering vector [p; vx; vy] for each VS
        idx = (n-1)*3 + (1:3);
        v_vs(idx) = [p_steer; vx_steer; vy_steer];
    end
end

% FUNCTION: PLOT BEAM PATTERN
% Generate beam pattern plot by scanning azimuthal angle around array centre
function [theta_deg, beam_pattern, fig_handle] = plot_beam_pattern(...
    response_db, X_grid, Y_grid, array_centre, plot_title, source_positions)

    % Reshape response to grid
    grid_response = reshape(response_db, size(X_grid));
    
    % Determine seach radius from source positions
    % Use median distance to sources
    distances = sqrt(sum((source_positions - array_centre').^2, 2));
    radius = median(distances);
    fprintf('  Using radius = %.3f m (median source distance)\n', radius);

    % Define angular resolution
    num_angles = 360; % 1 degree resolution
    theta_deg = linspace(0, 360, num_angles + 1);
    theta_deg = theta_deg(1:end-1); % Remove duplicate at 360
    theta_rad = deg2rad(theta_deg);
    
    % Calculate query points along each radial direction
    x_query = array_centre(1) + radius * sin(theta_rad);  
    y_query = array_centre(2) + radius * cos(theta_rad); %0deg is +Y direction
    
    % Interpolate beam response at query points
    beam_pattern = interp2(X_grid, Y_grid, grid_response, x_query, y_query, 'linear');
    
    % Handle any NaN values (points outside grid)
    if any(isnan(beam_pattern))
        warning('Some angles fall outside the scanned grid. Consider reducing radius or expanding scan area.');
        % Fill NaN with minimum value
        beam_pattern(isnan(beam_pattern)) = min(beam_pattern(~isnan(beam_pattern)));
    end
    
    % Normalise to 0 dB maximum
    beam_pattern = beam_pattern - max(beam_pattern);
    % Create figure with polar plot only
    fig_handle = figure('Name', plot_title, 'Position', [100, 100, 700, 700]);
    
    % Polar plot
    polarplot(theta_rad, beam_pattern, 'b-', 'LineWidth', 2);
    ax = gca;
    ax.ThetaDir = 'clockwise';
    ax.ThetaZeroLocation = 'top';
    rlim([min(beam_pattern), 0]);
    
    % Add radial grid labels
    rticks_vals = linspace(min(beam_pattern), 0, 5);
    rticks(rticks_vals);
    rticklabels(arrayfun(@(x) sprintf('%.1f dB', x), rticks_vals, 'UniformOutput', false));
    
    % Add source positions if provided
    if ~isempty(source_positions)
        hold on;
        for src = 1:size(source_positions, 1)
            % Calculate angle to source from array centre
            source_vec = source_positions(src, :)' - array_centre;
            source_angle_rad = atan2(source_vec(1), source_vec(2));
            % 0deg sits at positive y axis
            % Polar plot uses standard math convention
            r_lim = rlim;
            polarplot([source_angle_rad, source_angle_rad], r_lim, '--r', 'LineWidth', 2, ...
                'DisplayName', sprintf('Source %d', src));
        end
        legend('Location', 'northoutside', 'Orientation', 'horizontal');
        hold off;
    end
    
    % Set title and formatting
    title(sprintf('%s (r = %.3f m)', plot_title, radius), ...
        'FontName', 'Times New Roman', 'FontSize', 16, 'FontWeight', 'bold');
    ax.FontName = 'Times New Roman';
    ax.FontSize = 12;
    
    % Display beam width statistics
    fprintf('\n<strong>Beam Pattern Statistics:</strong>\n');
    fprintf('  Peak response: %.2f dB (at %.1f degrees)\n', 0, theta_deg(beam_pattern == max(beam_pattern)));
    
    % Calculate -3 dB beamwidth
    threshold_3db = max(beam_pattern) - 3;
    above_threshold = beam_pattern >= threshold_3db;
    
    % Find contiguous regions above threshold (ensure column vector)
    above_threshold = above_threshold(:); % Force to column vector
    diff_threshold = diff([0; above_threshold; 0]);
    starts = find(diff_threshold == 1);
    ends = find(diff_threshold == -1) - 1;
    
    if ~isempty(starts)
        fprintf('  -3 dB Beamwidth(s):\n');
        for i = 1:length(starts)
            if ends(i) >= starts(i)
                beamwidth = theta_deg(ends(i)) - theta_deg(starts(i));
                if beamwidth < 0
                    beamwidth = beamwidth + 360;
                end
                fprintf('    Region %d: %.1f degrees (%.1f° to %.1f°)\n', ...
                    i, beamwidth, theta_deg(starts(i)), theta_deg(ends(i)));
            end
        end
    end
    
    % Calculate sidelobe level
    main_lobe_mask = beam_pattern >= threshold_3db;
    main_lobe_mask = main_lobe_mask(:); % Force to column vector
    if any(~main_lobe_mask)
        beam_pattern_col = beam_pattern(:); % Also force beam_pattern to column
        sidelobe_level = max(beam_pattern_col(~main_lobe_mask));
        fprintf('  Maximum sidelobe level: %.2f dB\n', sidelobe_level);
    end
    
    % Calculate directivity index 
    % DI = 10*log10(4*pi / integral(pattern^2 dOmega))
    beam_linear = 10.^(beam_pattern(:)/10); % Force to column vector
    theta_rad_col = theta_rad(:); % Force to column vector
    integral_val = trapz(theta_rad_col, beam_linear);
    DI = 10*log10(2*pi / integral_val);
    fprintf('  Directivity Index: %.2f dB\n', DI);
    
    fprintf('\n');
end