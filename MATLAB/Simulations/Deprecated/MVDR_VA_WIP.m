clc; clear all; close all;

% L Marshall 11/11 
% Variable Number of Sources Localisation via an Acoustic Vector Array

% MVDR method of beamforming %
% Using broadband analysis with spectrum blocks %
% with eigenvalue decomposition, adaptive regularisation, and weighted MVDR

%% DEFINE INPUT VARIABLES %%

csv_filename = 'generatedsignal_avs.csv';

% Vector Sensor Array Parameters;
N_v = 2; %No. of vector sensors
d_y = 0.07; %Vector sensor y axis spacing (m)
delta = 0.0425; %MEMS colocation spacing (m)

% Define geom centre of the array 
x_a = 2.0 + delta/2; % Geometric centre x-coordinate
z_a = 0; % Array centre z-coordinate
y_a = -d_y/2 + delta/2; % Geometric centre y-coordinate

% Physical Environment Parameters
c_0 = 340; %Speed of sound (m/s)
rho_0 = 1.02; %Air Density at STP

% Sampling Characteristics
freq_limit = 3000; %Max freq for analysis (Hz)
overlap = 0.5; %Percentage overlap for sample window (%)
loading = 1e-4; %Set regularisation parameter (loading var in 2020 script)

% Grid Search Parameters
x_scan_points = 200; 
y_scan_points = 200;
x_margin = 3; %(m)
y_margin = 3; %(m)

% Visualisation Parameters
% Source positions [x, y, z] (m)
source_positions = [
    0, -0.25, 0; %Source 1
    %4, -2, 0.0; %Source 2
];

% Source frequencies (Hz)
source_frequencies = [
    1000; %Source 1 freq
    % 2000 %Source 2 freq
];

% Source amplituds 
source_amplitudes = [
    0.01; %Source 1 a
    % 0.01; %Source 2 a
];

num_sources = size(source_positions, 1);

%% Load in Generated Signal from CSV %%
% Data generated in 'VA_Output_WIP.m'

fprintf('Loading signal data from: %s\n', csv_filename);

if ~exist(csv_filename, 'file')
    error('Signal file not found: %s\nRun signal generation first.', csv_filename);
end

data_table = readtable(csv_filename);
data = table2array(data_table);
time = data(:,1); %Time vector
[numRows, numCols] = size(data);

N = numRows; %Number of samples
N_m = (numCols - 1) / 2; %Number of receiving elements 

fprintf('Loaded %d samples from %d microphones\n', N, N_m);

%% COMBINE SIGNAL COMPONENTS FOR PROCESSING %%
% Combine monopole and sinusoidal signals for each mic
tx_vs = zeros(N, N_m); %Preallocate transmitted signal array
for vec_sens = 1:N_v
    for mic = 1:4
        idx = (vec_sens-1)*4+mic;
        % Monopole Signals
        monopole_col = idx+1;
        % Sinusoidal Signals
        sinusoidal_col = 1+N_m+idx;

        tx_vs(:, idx) = data(:, monopole_col) + data(:, sinusoidal_col);
    end
end

tx_vs = tx_vs.'; %Transpose signals, rows = sensors, columns = time samples

%% ARRAY GEOM SET UP %%
% Adjust array geom to be centred 
if N_v == 1
    y_a_adjusted = y_a;
else
    y_a_adjusted = y_a - (N_v-1)*d_y/2;
end

% Define vector sensor colocation centres
vs_centres = zeros(3, N_v);
for vs = 1:N_v
    vs_centres(:, vs) = [x_a; y_a_adjusted + (vs-1)*d_y; z_a];
end

% Define all microphone positions 
mic_positions = zeros(3, N_m);
for vs = 1:N_v
    idx = (vs-1)*4 + (1:4);
    y_c = y_a_adjusted + (vs-1)*d_y;
    
    mic_positions(:,idx(1)) = [x_a - delta/2; y_c - delta/2; z_a]; %p0 - bottom left
    mic_positions(:,idx(2)) = [x_a + delta/2; y_c - delta/2; z_a]; %p1 - bottom right
    mic_positions(:,idx(3)) = [x_a + delta/2; y_c + delta/2; z_a]; %p2 - top right
    mic_positions(:,idx(4)) = [x_a - delta/2; y_c + delta/2; z_a]; %p3 - top left
end

%% CONVERT TO FREQ DOMAIN %%

fprintf('Converting to frequency domain...\n');

d_t = time(2) - time(1); %Time step (s)
F_s = 1/d_t; %Sampling frequency (Hz)
f = F_s * (0:floor(N/2)) / N; %Frequency vector

% Limit frequency range
idx_limit = f <= freq_limit;

% Preallocate FFT matrix
tx_freq = zeros(N_m, sum(idx_limit));

% Plot frequency spectrum
figure('Name', 'Frequency Domain');
hold on;
for sensor = 1:N_m
    % FFT of each sensor signal
    Y = fft(tx_vs(sensor, :));
    Y = Y(1:floor(N/2)+1); %Keep positive freqs (complex)
    % Store only the useful freqs (<= freq_limit)
    tx_freq(sensor,:) = Y(idx_limit);
    
    % Plot magnitude spectrum
    Y_mag = abs(Y)/N;
    Y_mag(2:end-1) = 2*Y_mag(2:end-1);
    Y_dB = 20*log10(Y_mag + eps);
    plot(f, Y_dB, 'DisplayName', sprintf('Sensor %d', sensor));
end
hold off;
xlabel('Frequency (Hz)', 'FontName', 'Times New Roman', 'FontSize', 14);
ylabel('Amplitude (dB)', 'FontName', 'Times New Roman', 'FontSize', 14);
xlim([0 freq_limit]);
legend('show', 'Location', 'best');
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);


%% DEFINING THE GRID SEARCH AREA %%

fprintf('Defining search area...\n');

x_scan = linspace(x_a - x_margin, x_a + x_margin, x_scan_points);
y_scan = linspace(y_a - y_margin, y_a + y_margin, y_scan_points);
z_scan = 0;  %2D plane change if needed

[X_grid, Y_grid] = meshgrid(x_scan, y_scan);
Z_grid = zeros(size(X_grid));

%N_mx3 matrix
candidate_points = [X_grid(:), Y_grid(:), Z_grid(:)];

% Grid resolution in x and y directions
x_res = x_scan(2) - x_scan(1);
y_res = y_scan(2) - y_scan(1);
grid_res = mean([x_res, y_res]);

fprintf('Search grid: %d x %d points, resolution: %.4f m\n', ...
    length(x_scan), length(y_scan), grid_res);

%% CREATING FREQ BINS FOR BROADBAND ANALYSIS %%

fprintf('Creating frequency bins...\n');

% Maximum bin width for req spatial resolution
max_binwidth = 1/(8*d_y*(N_m-1)/c_0); %(Hz)
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

snapshots_vs = make_snapshots(tx_vs, size_fft, overlap, window);
fprintf('  Created %d snapshots\n', size(snapshots_vs, 3));

r_vs = create_vs_csm(snapshots_vs, bin_index, delta, bin_freqs, rho_0, N_v, num_bins, c_0);
fprintf('  Created %dx%dx%d cross-spectral matrix\n', size(r_vs));

%% PERFORM BEAMFORMING %%

fprintf('\n<strong>Running MVDR Beamforming...</strong>\n');
fprintf('  Regularisation parameter: %.2e\n', loading);

% 1D scan along y-axis along set x value (for verification)
for src=1:num_sources
    fprintf('\n1D Scan (x = %.4f):\n', source_positions(src, 1));
    [scan_1d_max_db, scan_1d_y] = plot_1dscan_mvdr(r_vs, vs_centres, ...
        source_positions(1, 2), bin_freqs, c_0, loading, rho_0, num_bins, N_v, source_positions(src, 1));
end

% First calculate the response_db for the entire grid
response_db = mvdr_vs_beamforming(r_vs, vs_centres, candidate_points, ...
    bin_freqs, c_0, rho_0, loading, num_bins, N_v);

% 2D scan of full search area
fprintf('\n2D Scan:\n');
[est_x, est_y, scan_2d_max_db, all_est_positions] = plot_2dscan_mvdr(r_vs, vs_centres, ...
    candidate_points, bin_freqs, c_0, rho_0, y_scan, x_scan, ...
    X_grid, Y_grid, loading, num_bins, N_v, num_sources, source_positions, response_db);


%% PERFORMANCE METRICS %%

fprintf('\n<strong>Performance Metrics:</strong>\n');

% Initialise storage for all source results
all_results = struct();
all_results.est_positions = zeros(num_sources, 2);
all_results.true_positions = zeros(num_sources, 2);
all_results.errors = zeros(num_sources, 2);
all_results.MSE = zeros(num_sources, 1);
all_results.MSE_percent = zeros(num_sources, 1);
all_results.radial_error = zeros(num_sources, 1);
all_results.angular_error = zeros(num_sources, 1);
all_results.peak_response = zeros(num_sources, 1);

% Calculate errors relative to each source 
for src=1:num_sources
    fprintf('<strong>SOURCE %d: (%.3f, %.3f) m @ %.0f Hz</strong>\n', ...
        src, source_positions(src,1), source_positions(src,2), source_frequencies(src));

    % Define search region around this source
    search_radius_x = 0.5;  %(m)
    search_radius_y = 0.5;  %(m)

    % Extract true position for this source
    true_x = source_positions(src, 1);
    true_y = source_positions(src, 2);

    % Find grid points within search region
    x_mask = (X_grid >= true_x - search_radius_x) & (X_grid <= true_x + search_radius_x);
    y_mask = (Y_grid >= true_y - search_radius_y) & (Y_grid <= true_y + search_radius_y);
    region_mask = x_mask & y_mask;

    % Get beamforming response in this region
    grid_response = reshape(response_db, length(y_scan), length(x_scan));
    region_response = grid_response;
    region_response(~region_mask) = -Inf; %Mask out other regions

    % Find peak in this region
    [max_val, max_idx] = max(region_response(:));
    [y_idx, x_idx] = ind2sub(size(grid_response), max_idx);
    est_x = X_grid(y_idx, x_idx);
    est_y = Y_grid(y_idx, x_idx);

    % Calculate errors
    est_pos = [est_x, est_y];
    true_pos = [true_x, true_y];
    errors = true_pos - est_pos;

    squared_dist = sum(errors.^2);
    MSE = squared_dist;
    MSE_percent = (MSE / (grid_res^2)) * 100;
    radial_error = sqrt(squared_dist);

    % Angular error
    y_sensor = vs_centres(2,1);
    true_angle = atan2(true_y - y_sensor, true_x - x_a);
    est_angle = atan2(est_y - y_sensor, est_x - x_a);
    angular_error_deg = rad2deg(abs(true_angle - est_angle));
    % Normalise to 180 deg
    if angular_error_deg > 180
        angular_error_deg = mod(angular_error_deg + 180, 360) - 180;
    end

    % Store results
    all_results.est_positions(src,:) = est_pos;
    all_results.true_positions(src,:) = true_pos;
    all_results.errors(src,:) = errors;
    all_results.MSE(src) = MSE;
    all_results.MSE_percent(src) = MSE_percent;
    all_results.radial_error(src) = radial_error;
    all_results.angular_error(src) = angular_error_deg;
    all_results.peak_response(src) = max_val;


    % Display results for this source
    fprintf('Estimated position: (%.4f, %.4f) m\n', est_x, est_y);
    fprintf('Position error: dx = %.4f m, dy = %.4f m\n', errors(1), errors(2));
    fprintf('Radial error: %.4f m\n', radial_error);
    fprintf('Mean Squared Error: %.4f m²\n', MSE);
    fprintf('MSE (normalized): %.3f%% of grid resolution\n', MSE_percent);
    fprintf('Angular Error: %.3f deg\n', angular_error_deg);
    fprintf('Peak Response: %.3f dB\n', max_val);
    fprintf('\n');

end

% SUMMARY OF RESULTS
fprintf('<strong>OVERALL SUMMARY (%d sources)</strong>\n', num_sources);
fprintf('Average MSE: %.4f m² (± %.4f)\n', ...
    mean(all_results.MSE), std(all_results.MSE));
fprintf('Average radial error: %.3f m (± %.3f)\n', ...
    mean(all_results.radial_error), std(all_results.radial_error));
fprintf('Average angular error: %.4f deg (± %.4f)\n', ...
    mean(all_results.angular_error), std(all_results.angular_error));
fprintf('Best MSE: %.4f m² (Source %d)\n', ...
    min(all_results.MSE), find(all_results.MSE == min(all_results.MSE), 1));
fprintf('Worst MSE: %.4f m² (Source %d)\n', ...
    max(all_results.MSE), find(all_results.MSE == max(all_results.MSE), 1));
fprintf('\n');

%% FUNCTION DEFINITION %%

% FUNCTION: MAKE SNAPSHOTS OF GENERATED SIGNAL
% Make time domain snapshots of input signal
function snapshots = make_snapshots(tx_vs, size_fft, overlap, window)
    % No. sample to move per snapshot
    step = round(overlap * size_fft);
    % No. of full snapshots that fit in signal with given step size
    num_snap = floor((size(tx_vs,2) - size_fft) / step) + 1; %No. snapshots
    % Start indices for each snapshot
    start_idx = 1 + (0:(num_snap-1)) * step;
    
    % Create index matrix (size_fft x num_snap)
    idx_matrix = start_idx + (0:size_fft-1)'; %each column = indices for a snapshot
    % Extract snapshots (force reshape into 3D)
    snapshots = tx_vs(:, idx_matrix);
    snapshots = reshape(snapshots, size(tx_vs,1), size_fft, num_snap);
    % Apply window along snapshots dimension
    snapshots = snapshots .* reshape(window, [1, size_fft, 1]);
end


% FUNCTION: FFT AND CREATE CROSS-SPECTRAL MATRIX OF VECTOR SENSORS
%Convert each block to freq domain and create cross-spectral matrix
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


% FUNCTION: MVDR BEAMFORMING
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

% FUNCTION: 2D SCAN PLOTTING
function [est_x, est_y, scan_2d_max_db, all_est_positions] = plot_2dscan_mvdr(...
    r_vs, vs_centres, candidate_points, bin_freqs, c_0, rho_0, y_scan, x_scan, ...
    X_grid, Y_grid, loading, num_bins, N_v, num_sources, source_positions, response_db)
    
    % If response_db not provided, calculate it
    if nargin < 18
        response_db = mvdr_vs_beamforming(r_vs, vs_centres, candidate_points, ...
            bin_freqs, c_0, rho_0, loading, num_bins, N_v);
    end
    
    % Reshape responses to a grid for plotting
    grid_response = reshape(response_db, length(y_scan), length(x_scan));
    
    figure('Name', '2D MVDR Beamforming - All Sources');
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
    est_x = X_grid(y_idx, x_idx);
    est_y = Y_grid(y_idx, x_idx);
    
    % Find estimated position for each source
    all_est_positions = zeros(num_sources, 2);
    search_radius_x = 0.5;
    search_radius_y = 0.5;
    
    source_colors = lines(num_sources);
    
    for src = 1:num_sources
        true_x = source_positions(src, 1);
        true_y = source_positions(src, 2);
        
        % Define search region
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
        plot(true_x, true_y, 'x', 'Color', source_colors(src,:), ...
            'MarkerSize', 14, 'LineWidth', 3, ...
            'DisplayName', sprintf('True Source %d', src));
        
        % Plot estimated position
        plot(est_x_src, est_y_src, 'o', 'Color', source_colors(src,:), ...
            'MarkerSize', 12, 'LineWidth', 2, 'MarkerFaceColor', 'none', ...
            'DisplayName', sprintf('Est Source %d', src));
        
    end
    
    legend('Location', 'northeast', 'Color', 'w', 'TextColor', 'k');
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
    title(sprintf('%d Sources: MVDR Beamforming', num_sources), ...
        'FontName', 'Times New Roman', 'FontSize', 14);
end


% FUNCTION: 1D SCAN PLOTTING
function [max_response, y_est] = plot_1dscan_mvdr(r_vs, vs_centres, ...
    source_y, bin_freqs, c_0, loading, rho_0, num_bins, N_v, x_fixed)
    
    y_scan_line = linspace(source_y-2, source_y+2, 200);
    %x_fixed = 0;
    z_fixed = 0;
    
    candidate_points_line = [x_fixed*ones(length(y_scan_line),1), ...
        y_scan_line(:), z_fixed*ones(length(y_scan_line),1)];
    
    response_line = mvdr_vs_beamforming(r_vs, vs_centres, ...
        candidate_points_line, bin_freqs, c_0, rho_0, loading, num_bins, N_v);
    
    [max_response, max_idx] = max(response_line);
    y_est = y_scan_line(max_idx);
    
    fprintf('  Max response: %.2f dB at y = %.3f m\n', max_response, y_est);
    
    figure('Name', '1D MVDR Scan');
    plot(y_scan_line, response_line, 'LineWidth', 2);
    xlabel('y position (m)', 'FontName', 'Times New Roman', 'FontSize', 14);
    ylabel('Beamformer Response (dB)', 'FontName', 'Times New Roman', 'FontSize', 14);
    grid on;
    hold on;
    
    plot(source_y, max(response_line), 'kx', 'MarkerSize', 10, 'LineWidth', 1.5);
    xline(source_y, 'r--', 'LineWidth', 1.5, 'DisplayName', 'True Source');
    
    legend('Beamformer Response', 'True Source');
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
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
    % Output: [p1; vx1; vy1; p2; vx2; vy2]
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
        % rho_c = (rho_0*c_0)*((1-1i/(k*r))/(1+1/(k*r)^2));
        rho_c = rho_0*c_0;
        vx_steer = rho_c * vx_steer;
        vy_steer = rho_c * vy_steer;
        % Store in steering vector [p; vx; vy] for each VS
        idx = (n-1)*3 + (1:3);
        v_vs(idx) = [p_steer; vx_steer; vy_steer];
    end
end

