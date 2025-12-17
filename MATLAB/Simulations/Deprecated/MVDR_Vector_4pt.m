clc; clear all; close all;

% L Marshall 06/11 
% Single Source Localisation via an Acoustic Vector Array

% MVDR method of beamforming %
% Using broadband analysis with spectrum blocks %
% with eigenvalue decomposition, adaptive regularisation, and weighted MVDR

%% Load in Generated Signal from CSV %%
% Data generated in 'Vector_Array_Output_0510.m'

data_table = readtable('generatedsignal_avs.csv');
data = table2array(data_table);
time = data(:,1); %Time vector
[numRows, numCols] = size(data);

N = numRows; %Number of samples
N_m = (numCols - 1)/2; %Number of receiving elements 
N_v = N_m/4; %Since 4 colocated components selected

% Combine monopole and sinusoidal signals
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

tx_vs = tx_vs.'; %Transpose signals

%% Initialising Array and Signal Parameters %%

% Inputs Inferred from Array Output Script
c_0 = 340; %Speed of sound (m/s)
rho_0 = 1.02; %Air density

d_y = 0.04; %Spacing between vector elements
delta = 0.04; %Finite diff spacing

x_a = 2.0 + delta/2; %Array geometric centre x (m)
y_a = -d_y/2 + delta/2; %Array geometric centre y (m)
z_a = 0; %Array centre z (m)


% Defining array centre positions
% vs_centres = [
%     x_a, -d_y/2, z_a;  % MEMS1 at y = -0.15
%     x_a,  d_y/2, z_a   % MEMS2 at y = +0.15
% ]';

vs_centres = zeros(3, N_v);
for vs = 1:N_v
  
    vs_centres(:, vs) = [x_a; y_a + (vs-1)*d_y; z_a];
end

%vs_centres = [x_a, -d_y/2, z_a];

% Reference source position
source_x = 0;
source_y = -0.25;

% Define all microphone positions (to match signal gen script)
mic_positions = zeros(3, N_m);
for vs = 1:N_v
    idx = (vs-1)*4 + (1:4);
    % y_c = vs_centres(2, vs);
    %y_c = -d_y/2;
    
    % Use x_a and y_a directly to match signal generation
    mic_positions(:,idx(1)) = [x_a - delta/2; y_a - delta/2; z_a];  % p0: bottom-left
    mic_positions(:,idx(2)) = [x_a + delta/2; y_a - delta/2; z_a];  % p1: bottom-right
    mic_positions(:,idx(3)) = [x_a + delta/2; y_a + delta/2; z_a];  % p2: top-right
    mic_positions(:,idx(4)) = [x_a - delta/2; y_a + delta/2; z_a];  % p3: top-left
end

%% Convert to Frequency Domain %%
% Fast fourier transform
d_t = time(2) - time(1); %Time step (s)
F_s = 1/d_t; %Sampling freq (Hz)
f = F_s*(0:floor(N/2))/N; %Freq vector

% Limit frequency range to 0–3000 Hz
freq_limit = 3000;
idx_limit = f <= freq_limit; %Index limit boolean check

% Preallocate matrix: rows = sensors, columns = frequency bins
tx_freq = zeros(N_m, sum(idx_limit));

figure;
hold on;
for sensor = 1:N_m
    % FFT of each sensor signal
    Y = fft(tx_vs(sensor, :));
    Y = Y(1:floor(N/2)+1); %Keep positive freqs (complex)
    
    % Store only the useful freqs (<= 1000 Hz)
    tx_freq(sensor,:) = Y(idx_limit);
    
    % Plot magnitude spectrum 
    Y_mag = abs(Y)/N;
    Y_mag(2:end-1) = 2*Y_mag(2:end-1);
    Y_dB = 20*log10(Y_mag);
    plot(f, Y_dB, 'DisplayName', sprintf('Sensor %d', sensor));
end
hold off;
xlabel('Frequency (Hz)');
ylabel('Amplitude (dB)');
% title('Transmitted Signal in Frequency Domain');
xlim([0 freq_limit]);
legend show;
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14); % optional for report styling

%% Defining the Search Area %%
% Note: need to make it around the array centre
x_margin = 3;
y_margin = 1;

x_scan = linspace(x_a-x_margin, x_a+x_margin, 100); 
y_scan = linspace(y_a-y_margin, y_a+y_margin, 100);

z_scan = 0; %2D plane for now

[X_grid, Y_grid] = meshgrid(x_scan, y_scan);
Z_grid = zeros(size(X_grid));

candidate_points = [X_grid(:), Y_grid(:), Z_grid(:)]; %N_mx3 matrix

% Grid resolution in x and y directions
x_res = x_scan(2) - x_scan(1);
y_res = y_scan(2) - y_scan(1);

grid_res = mean([x_res, y_res]);


%% Creating Frequency Bins for Broadband Analysis %%

% Max bin width to achieve req. spatial res
max_binwidth = 1/(8*d_y*(N_m-1)/c_0); %(Hz)
size_fft = floor(F_s/max_binwidth); %No. points that will be processed in signal

% Sampling frequency and FFT size
delta_f = F_s/N; %Frequency spacing per bin

% Define frequency vector
fft_vec = F_s*(0:(size_fft-1))/size_fft; 

% Target frequencies based on required spatial res
% Step in terms of FFT bins
min_freq = max_binwidth;
target_freqs = min_freq:max_binwidth:freq_limit;

% Map each target to nearest snapshot FFT bin
bin_index = zeros(size(target_freqs));
for k = 1:length(target_freqs)
    [~,bin_index(k)] = min(abs(fft_vec-target_freqs(k)));
end

% Generate bin indices within accepted range 0–3000 Hz
bin_index = unique(bin_index); %Remove duplicates
bin_index = bin_index(bin_index <= size_fft);
bin_freqs = fft_vec(bin_index); %Actual bin freqs
num_bins = length(bin_index); %No. bins

% Define sampling characteristics
overlap = 0.5; %Percentage overlap (50%)
num_snaps = 8*N_m; %No. snapshots

% Calc no. samples accounting for overlap
num_samps = (num_snaps+1)*size_fft*overlap;

% Define a taper for sampling window
window = hanning(size_fft)';

%% Create AVS snapshots and CSM %%
snapshots_vs = make_snapshots(tx_vs, size_fft, overlap, window);
fprintf('Created %d snapshots from signal\n', size(snapshots_vs, 3));
r_vs = create_vs_csm(snapshots_vs, bin_index, delta, bin_freqs, rho_0, N_v, num_bins, c_0);

%% Perform Beamforming %%

fprintf('<strong>Acoustic Vector Sensor MVDR Beamforming in Broadband Spectrum</strong>\n')

% Set regularisation parameter (aka loading var from 2020 script)
loading = 1e-4; 

% Run MVDR beamforming
%response_db = mvdr_beamforming(r, mic_positions, candidate_points, bin_freqs, c_0, loading);
response_db = mvdr_vs_beamforming(r_vs, vs_centres, candidate_points, ...
    bin_freqs, c_0, rho_0, loading, num_bins, N_v);

% Plot the 1D scan of beamformed response 
% along x=0 to check correct peak found
[scan_1d_max_db, scan_1d_y] = plot_1dscan_mvdr(r_vs, vs_centres, source_y, ...
    bin_freqs, c_0, loading, rho_0, num_bins, N_v);

% Plot the 2D scan of beamformed response
[est_x, est_y, max_val] = plot_2dscan_mvdr(r_vs, vs_centres, candidate_points, ...
    bin_freqs, c_0, rho_0, y_scan, x_scan, source_x, source_y, X_grid, ...
    Y_grid, loading, num_bins, N_v);

% Display performance metrics
est_pos = [est_x est_y];
true_pos = [source_x source_y];
errors = true_pos - est_pos;
squared_dist = sum(errors.^2, 2); %Squared distance per source
MSE = mean(squared_dist); %Average across sources
% Normalised MSE w.r.t grid resolution
MSE_percent = (MSE / (grid_res^2)) * 100;
y_sensor = vs_centres(2,1);
% True DOA angle
true_angle = atan2(source_y - y_sensor, source_x - x_a);
% Estimated DOA angle  
est_angle = atan2(est_y - y_sensor, est_x - x_a);
% Angular error
angular_error_deg = rad2deg(abs(true_angle - est_angle));

fprintf('\n<strong>Performance Outputs:</strong>\n');
fprintf('Mean Squared Error: %.4f m^2\n', MSE);
fprintf('Mean Squared Error (normalised): %.2f %% of avg. grid resolution\n', MSE_percent);
fprintf('Estimation Error to True Position: dx = %.3f m, dy = %.3f m\n', errors(1), errors(2));
fprintf('Angular Error: %.4f deg\n', angular_error_deg);
%fprintf('Radial Error Distance: %.3f m\n', sqrt(sum(errors.^2)));

%% Function Definition %%

% FUNCTION ONE: MAKE SNAPSHOTS OF GENERATED SIGNAL
% Make time domain snapshots of input signal
function snapshots = make_snapshots(tx_vs, size_fft, overlap, window)
    % No. sample to move per snapshot
    step = round(overlap*size_fft);
    
    % No. of full snapshots that fit in signal with given step size
    num_snap = floor((size(tx_vs,2)-size_fft)/step)+1; %No. snapshots 
    % Start indices for each snapshot
    start_idx = 1+(0:(num_snap-1))*step;
    
    % Create index matrix (size_fft x num_snap)
    idx_matrix = start_idx+(0:size_fft-1)';  %each column = indices for a snapshot
    
    % Extract snapshots (force reshape into 3D)
    snapshots = tx_vs(:, idx_matrix);  
    snapshots = reshape(snapshots, size(tx_vs,1), size_fft, num_snap);
    
    % Apply window along snapshots dimension
    snapshots = snapshots.*reshape(window, [1, size_fft, 1]);
end


% FUNCTION TWO: FFT AND CREATE CROSS-SPECTRAL MATRIX OF VECTOR SENSORS%
% Convert each block to freq domain and create cross-spectral matrix
function r_vs = create_vs_csm(snapshots_vs, bin_index, delta, bin_freqs, rho_0, N_v, num_bins, c_0)
    % Create CSM for AVS array
    % snapshots: 6 x size_fft x num_snap (6 mics total)
    
    %num_bins = length(bin_index);
    signal_dim = 3*N_v; %(p, vx, vy for each AVS)
    
    % FFT of snapshots
    fx = fft(snapshots_vs, [], 2);
    fx = fx(:, bin_index, :);
    
    % Initialize CSM
    r_vs = zeros(signal_dim, signal_dim, num_bins);

    num_snap = size(snapshots_vs, 3);
    
    for jf = 1:num_bins
        freq = bin_freqs(jf);
        
        % Average over snapshots
        fx_bin = squeeze(fx(:, jf, :)); %M_mics x N_snap
        %fx_avg = mean(fx_bin, 2);

        % Initialise CSM for this frequency bin
        R_freq = zeros(signal_dim, signal_dim);

        % Loop over snapshots
        for snap = 1:num_snap
            fx_snap = fx_bin(:, snap); % Pressure at each mic for this snapshot
            
            % Compute VS outputs for THIS snapshot
            [p_vs, vx_vs, vy_vs] = vs_outputs(fx_snap, delta, freq, rho_0, N_v, c_0);
            
            % Apply impedance scaling (plane wave approximation)
            rho_c = rho_0*c_0;
            vx_vs = rho_c*vx_vs;
            vy_vs = rho_c*vy_vs;
            
            % Form signal vector for this snapshot
            vs_signal = [p_vs(1); vx_vs(1); vy_vs(1)]; 
            
            % Accumulate outer product
            R_freq = R_freq + (vs_signal * vs_signal');
        end
        
        % Average over snapshots
        r_vs(:,:,jf) = R_freq / num_snap;
    end
end

% FUNCTION THREE: MVDR BEAMFORMING %
% MVDR beamforming using methods from Dr.Bao script
function response_db = mvdr_vs_beamforming(r_vs, vs_centres, candidate_points, ...
    bin_freqs, c_0, rho_0, loading, num_bins, N_v)
    
    num_points = size(candidate_points, 1);
    %num_bins = size(r_vs, 3);
    
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

% FUNCTION FOUR: 2D SCAN PLOTTING FOR MVDR %
function [est_x, est_y, max_val] = plot_2dscan_mvdr(r_vs, vs_centres, candidate_points, ...
    bin_freqs, c_0, rho_0, y_scan, x_scan, source_x, source_y, ...
    X_grid, Y_grid, loading, num_bins, N_v)
    
    fprintf('Creating 2D plot.\n');

    % Run advanced MVDR beamforming
    response_db = mvdr_vs_beamforming(r_vs, vs_centres, candidate_points, ...
    bin_freqs, c_0, rho_0, loading, num_bins, N_v);

    % Reshape responses to a grid for plotting
    grid_response = reshape(response_db, length(y_scan), length(x_scan));

    % Plot spatial response map
    figure;
    imagesc(x_scan, y_scan, grid_response);
    axis xy;
    xlabel('X (m)');
    ylabel('Y (m)');
    % title('MVDR Beamforming Response');
    colorbar;
    colormap('jet');
    hold on;

    % Find estimated source coordinates
    [max_val, max_idx] = max(grid_response(:));
    [y_idx, x_idx] = ind2sub(size(grid_response), max_idx);
    est_x = X_grid(y_idx, x_idx);
    est_y = Y_grid(y_idx, x_idx);

    % For debugging
    fprintf('2D MVDR max: %.2f dB at x = %.3f m, y = %.3f m\n', max_val, est_x, est_y);

    % Plot true source (red x) and estimated source (blue circle)
    plot(source_x, source_y, 'rx', 'MarkerSize', 12, 'LineWidth', 2);
    plot(est_x, est_y, 'bo', 'MarkerSize', 10, 'LineWidth', 2);
    
    legend({'True Source', 'Estimated'}, ...
        'Location', 'northeast', 'Color', 'w', 'TextColor', 'k');
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);

end

% FUNCTION FIVE: 1D SCAN PLOTTING FOR MVDR
function [max_response, y_est] = plot_1dscan_mvdr(r_vs, vs_centres, source_y, bin_freqs, c_0, loading, ...
    rho_0, num_bins, N_v)

    fprintf('Creating 1D plot for peak testing.\n');
    % Scan along y-axis at fixed x and z
    y_scan_line = linspace(-0.5, 0.5, 200);
    x_fixed = 0; 
    z_fixed = 0;

    % Build candidate points
    candidate_points_line = [x_fixed*ones(length(y_scan_line),1), ...
        y_scan_line(:), z_fixed*ones(length(y_scan_line),1)];

    % Run advanced MVDR beamforming along the line
    response_line = mvdr_vs_beamforming(r_vs, vs_centres, ...
        candidate_points_line, bin_freqs, c_0, rho_0, loading, num_bins, N_v);

    % Find maximum response
    [max_response, max_idx] = max(response_line);
    y_est = y_scan_line(max_idx);

    % For debugging
    fprintf('1D Scan at x=0: max response = %.2f dB at y = %.3f m\n', max_response, y_est);
    
    % Plot the line response
    figure;
    plot(y_scan_line, response_line, 'LineWidth', 2);
    xlabel('y position (m)');
    ylabel('Beamformer Response (dB)');
    % title('1D Scan of MVDR Beamforming along y-axis');
    grid on;
    hold on;

    % Plot true source and estimated source
    plot(source_y, max(response_line), 'kx', 'MarkerSize', 10, 'LineWidth', 1.5); 
    xline(source_y, 'r--', 'LineWidth', 1.5, 'DisplayName','True Source');

    legend('Beamformer Response','True Source');
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 14); 
end


% % FUNCTION SIX: COMPUTE VECTOR SENSOR PRESSURE - VELOCITY CONVERSIONS
% function [p_vs, vx_vs, vy_vs] = vs_outputs(tx_freq, delta, freq, rho_0, N_v, c_0)
%     omega = 2*pi*freq;
%     k = 2*pi*freq/c_0;
%     % Initialise outputs
%     p_vs = zeros(N_v, 1);
%     vx_vs = zeros(N_v, 1);
%     vy_vs = zeros(N_v, 1);
% 
%     for vs = 1:N_v
%         idx = (vs-1)*4+(1:4);
% 
%         % Get pressure at 4 microphone positions (square corners)
%         p0 = tx_freq(idx(1)); %bottom-left: (x_a, y_c)
%         p1 = tx_freq(idx(2)); %bottom-right: (x_a+delta, y_c)
%         p2 = tx_freq(idx(3)); %top-right: (x_a+delta, y_c+delta)
%         p3 = tx_freq(idx(4)); %top-left: (x_a, y_c+delta)
% 
%         % Pressure estimate (average of all 4 corners)
%         p_vs(vs) = (p0 + p1 + p2 + p3)/4;
% 
%         % Particle velocity using 4-point square finite difference
%         % v = -delta_p/(i*omega*rho_0)
% 
%         % vx: average of bottom and top x-differences
%         % dp/dx ≈ [(p1-p0) + (p2-p3)]/(2delta)
%         vx_vs(vs) = (p0 + p3 - p1 - p2)/(2*1i*omega*rho_0*delta);
% 
%         % vy: average of left and right y-differences  
%         % dp/dy ≈ [(p3-p0) + (p2-p1)]/(2delta)
%         vy_vs(vs) = (p0 + p1 - p3 - p2)/(2*1i*omega*rho_0*delta);
%     end
% end

% FUNCTION SIX: COMPUTE VECTOR SENSOR OUTPUTS USING LEAST-SQUARES GRADIENT
function [p_vs, vx_vs, vy_vs] = vs_outputs(tx_freq, delta, freq, rho_0, N_v, c_0)
    omega = 2*pi*freq;
    
    % Initialise outputs
    p_vs = zeros(N_v, 1);
    vx_vs = zeros(N_v, 1);
    vy_vs = zeros(N_v, 1);

    for vs = 1:N_v
        idx = (vs-1)*4+(1:4);

        % Get pressure at 4 microphone positions (square corners)
        p0 = tx_freq(idx(1)); %bottom-left
        p1 = tx_freq(idx(2)); %bottom-right
        p2 = tx_freq(idx(3)); %top-right
        p3 = tx_freq(idx(4)); %top-left
        
        % Microphone positions relative to centre
        % Centre is at (0, 0) in local coords
        mic_positions = [
            -delta/2, -delta/2;  %p0: bottom-left
             delta/2, -delta/2;  %p1: bottom-right
             delta/2,  delta/2;  %p2: top-right
            -delta/2,  delta/2   %p3: top-left
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

% FUNCTION SEVEN: STEERING VECTOR FOR ACOUSTIC VECTOR ARRAY
function v_vs = vs_steering_vector(vs_centres, source_pos, freq, c_0, rho_0, N_v)
    % Output: [p1; vx1; vy1; p2; vx2; vy2]
    
    k = 2*pi*freq/c_0;
    omega = 2*pi*freq;
    
    v_vs = zeros(3*N_v, 1); %3 components per vector sensor
    
    for n = 1:N_v
        % Vector from source to VS centre
        r_vec = vs_centres(:,n) - source_pos;
        r = norm(r_vec);
        r_hat = r_vec / r; %Unit vector
        
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

