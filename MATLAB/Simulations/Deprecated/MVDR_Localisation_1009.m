clc; clear all; close all;

% L Marshall 15/09 
% Single Source Localisation via a Linear Microphone Array

% MVDR method of beamforming %
% Using broadband analysis with spectrum blocks %
% with eigenvalue decomposition, adaptive regularisation, and weighted MVDR

%% Load in Generated Signal from CSV %%
% Data generated in 'Array_Output_0409.m'

data = readmatrix('generatedsignal.csv');
time = data(:,1); %Time vector
[numRows, numCols] = size(data);

N = numRows; %Number of samples
Nr = (numCols - 1)/2; %Number of receiving elements 

% Preallocate transmitted signal array
tx = zeros(N, Nr); 

for mic = 1:Nr
    % Monopole Signals
    monopole_col = mic + 1; %Columns 2-12
    % Sinusoidal Signals
    sinusoidal_col = mic + 12; %Columns 13-23

    tx(:, mic) = data(:, monopole_col) + data(:, sinusoidal_col);
end

tx = tx.'; %Transpose signals

%% Initialising Array and Signal Parameters %%

% Inputs Inferred from Array Output Script
c_0 = 340; %Speed of sound (m/s)
d_y = 0.08; %Sensor spacing (m)
x_a = 2.0; %Array centre x coord (m)
z_a = 0; %Array centre z coord (m)
y_a = 0; %Array centre y coord (m)

% For reference
source_x = 0;
source_y = -0.25;

n = linspace(-(Nr-1)/2, (Nr-1)/2, Nr); %No. of sensors
mic_positions = [x_a*ones(1,Nr); d_y*n; zeros(1,Nr)]; %Sensor positions
 
%% Convert to Frequency Domain %%
% Fast fourier transform
d_t = time(2) - time(1); %Time step (s)
F_s = 1/d_t; %Sampling freq (Hz)
f = F_s*(0:floor(N/2))/N; %Freq vector

% Limit frequency range to 0–3000 Hz
freq_limit = 3000;
idx_limit = f <= freq_limit; %Index limit boolean check

% Preallocate matrix: rows = sensors, columns = frequency bins
tx_freq = zeros(Nr, sum(idx_limit));

figure;
hold on;
for sensor = 1:Nr
    % FFT of each sensor signal
    Y = fft(tx(sensor, :));
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
% Note: need to make it around the array center
x_margin = 3;
y_margin = 1;

x_scan = linspace(x_a-x_margin, x_a+x_margin, 100); 
y_scan = linspace(y_a-y_margin, y_a+y_margin, 100);

z_scan = 0; %2D plane for now

[X_grid, Y_grid] = meshgrid(x_scan, y_scan);
Z_grid = zeros(size(X_grid));

candidate_points = [X_grid(:), Y_grid(:), Z_grid(:)]; %Nrx3 matrix

% Grid resolution in x and y directions
x_res = x_scan(2) - x_scan(1);
y_res = y_scan(2) - y_scan(1);

grid_res = mean([x_res, y_res]);


%% Creating Frequency Bins for Broadband Analysis %%

% Max bin width to achieve req. spatial res
max_binwidth = 1/(8*d_y*(Nr-1)/c_0); %(Hz)

% Time domain res delta_f = 1/block_length
% Calc min overall block length accounting for spatial res
%maxlen_block = 1/max_binwidth; %(secs)
size_fft = floor(F_s/max_binwidth); %No. points that will be processed in signal

% Sampling frequency and FFT size
delta_f = F_s/N; %Frequency spacing per bin

% Define frequency vector
fft_vec = F_s*(0:(size_fft-1))/size_fft; 

% Target frequencies based on required spatial res
% Step in terms of FFT bins
target_freqs = 0:max_binwidth:freq_limit;

% Map each target to nearest snapshot FFT bin
bin_index = zeros(size(target_freqs));
for k = 1:length(target_freqs)
    [~,bin_index(k)] = min(abs(fft_vec-target_freqs(k)));
end

% Generate bin indices within accepted range 0–1000 Hz
bin_index = unique(bin_index); %Remove duplicates
bin_index = bin_index(bin_index <= size_fft);
bin_freqs = fft_vec(bin_index); %Actual bin freqs
num_bins = length(bin_index); %No. bins

% Define sampling characteristics
overlap = 0.5; %Percentage overlap (50%)
num_snaps = 4*Nr; %No. snapshots

% Calc no. samples accounting for overlap
num_samps = (num_snaps+1)*size_fft*overlap;

% Define a taper for sampling window
window = hanning(size_fft)';

% Call on snapshot making function
snapshots = make_snapshots(tx, size_fft, overlap, window);

% Max no. snapshots that signal fits
maxnum_snap = size(snapshots,3); 

% Call function to fft and create cross-spectral matrix avg. over snapshots
r = create_csm(snapshots, bin_index, Nr, num_bins, maxnum_snap);


%% Perform Beamforming %%

fprintf('<strong>MVDR Beamforming in Broadband Spectrum</strong>\n')

% Set regularisation parameter (aka loading var from 2020 script)
loading = 1e-4; 

% Run MVDR beamforming
response_db = mvdr_beamforming(r, mic_positions, candidate_points, bin_freqs, c_0, loading);

% Plot the 1D scan of beamformed response 
% along x=0 to check correct peak found
plot_1dscan_mvdr(r, mic_positions, source_y, bin_freqs, c_0, loading)

% Plot the 2D scan of beamformed response
[est_x, est_y] = plot_2dscan_mvdr(r, mic_positions, candidate_points, bin_freqs, c_0, ...
    y_scan, x_scan, source_x, source_y, X_grid, Y_grid, loading);

% Display performance metrics
est_pos = [est_x est_y];
true_pos = [source_x source_y];
errors = true_pos - est_pos;
squared_dist = sum(errors.^2, 2); %Squared distance per source
MSE = mean(squared_dist); %Average across sources
% Normalised MSE w.r.t grid resolution
MSE_percent = (MSE / (grid_res^2)) * 100;

fprintf('\n<strong>Performance Outputs:</strong>\n');
fprintf('Mean Squared Error: %.4f m^2\n', MSE);
fprintf('Mean Squared Error (normalised): %.2f %% of avg. grid resolution\n', MSE_percent);
fprintf('Estimation Error to True Position: dx = %.3f m, dy = %.3f m\n', errors(1), errors(2));
%fprintf('Radial Error Distance: %.3f m\n', sqrt(sum(errors.^2)));

%% Function Definition %%

% FUNCTION ONE: MAKE SNAPSHOTS OF GENERATED SIGNAL
% Make time domain snapshots of input signal
function snapshots = make_snapshots(tx, size_fft, overlap, window)
    % No. sample to move per snapshot
    step = round(overlap*size_fft);
    
    % No. of full snapshots that fit in signal with given step size
    num_snap = floor((size(tx,2)-size_fft)/step)+1; %No. snapshots 
    % Start indices for each snapshot
    start_idx = 1+(0:(num_snap-1))*step;
    
    % Create index matrix (size_fft x num_snap)
    idx_matrix = start_idx+(0:size_fft-1)';  %each column = indices for a snapshot
    
    % Extract snapshots (force reshape into 3D)
    snapshots = tx(:, idx_matrix);  
    snapshots = reshape(snapshots, size(tx,1), size_fft, num_snap);
    
    % Apply window along snapshots dimension
    snapshots = snapshots.*reshape(window, [1, size_fft, 1]);
end


% FUNCTION TWO: FFT AND CREATE CROSS-SPECTRAL MATRIX %
% Convert each block to freq domain and create cross-spectral matrix
function r = create_csm(snapshots, bin_index, Nr, num_bins, maxnum_snap)
    % FFT along the time dimension for all snapshots
    fx = fft(snapshots,[],2); %(Nr x size_fft x num_snap)    
    
    % Keep only the frequencies of interest
    fx = fx(:,bin_index,:); %(Nr x num_bins x num_snap)
    
    % Preallocate CSM
    r = zeros(Nr,Nr,num_bins);
    
    % Vectorised additions to CSM per snapshot
    for jf = 1:num_bins
        % Extract matrix for this frequency bin
        fx1 = squeeze(fx(:,jf,:)); %(Nrxnum_snap)
    
        % Compute CSM: sum over snapshots of outer products
        r(:,:,jf) = fx1*fx1'; %(NrxNr)
    end
    
    % Average CSM over snapshots of sampled signal
    r = r/maxnum_snap;
end


% FUNCTION THREE: MVDR BEAMFORMING %
% MVDR beamforming using methods from Dr.Bao script
function response_db = mvdr_beamforming(r, mic_positions, candidate_points, bin_freqs, c_0, loading)
% 
% INPUTS:
% r: Nr x Nr x num_bins cross-spectral matrix
% mic_positions: 3 x Nr sensor positions
% candidate_points: N_points x 3 scan points

num_points = size(candidate_points, 1); %No. points to search
num_bins = size(r, 3);
Nr = size(r, 1);

% Preallocate arrays
mvdr_responses = zeros(num_points, num_bins);
cbf_responses = zeros(num_points, num_bins);


for jf = 1:num_bins
    k = 2*pi*bin_freqs(jf)/c_0; %Wave no. for this freq bin
    rf = squeeze(r(:,:,jf)); %Covariance matrix for this freq bin
    
    % Eigenvalue decomposition (improve stability)
    [ur, er] = eig(rf);
    erv = real(diag(er));
    
    % Safety check for negative eigenvalues 
    erv = max(erv, 1e-12);
    
    for n = 1:num_points
        l_pt = candidate_points(n,:).';
        
        % Calculate steering vector
        r_lm = sqrt(sum((mic_positions - l_pt).^2, 1));
        v1 = exp(-1i*k*r_lm).'; %Steering vector 
        
        % Conventional DAS beamformer output 
        cbf_out = real(v1' * rf * v1);
        cbf_responses(n, jf) = cbf_out;
        
        % Adaptive regularisation based on CBF output
        lambda = loading * cbf_out; %lambda = diagonal loading factor
        if lambda == 0
            lambda = loading; %Default to fixed loading
        end
        
        % Weighted MVDR using eigendecomposition
        rxv = (ur*diag(1./(erv+lambda))*ur')*v1;
        
        % Normalise to get weight vector
        denominator = v1'*rxv;
        if abs(denominator) < 1e-12
            mvdr_responses(n, jf) = 0;
        else
            vmvdr = rxv/denominator; %MVDR weight vector
            
            % Calc weighted output power
            mvdr_responses(n, jf) = abs(vmvdr'*rf*vmvdr);
        end
    end
    
    % Progress indicator
    % if mod(jf, max(1, floor(num_bins/10))) == 0
    %     fprintf('Processed %d/%d frequency bins\n', jf, num_bins);
    % end
end

% Broadband response by summing across frequencies
responses_sum = sum(mvdr_responses, 2);

% Convert to dB and normalise
response_db = 10*log10(responses_sum + eps);
response_db = response_db - max(response_db);
end


% FUNCTION FOUR: 2D SCAN PLOTTING FOR MVDR %
function [est_x, est_y] = plot_2dscan_mvdr(r, mic_positions, candidate_points, ...
    bin_freqs, c_0, y_scan, x_scan, source_x, source_y, X_grid, Y_grid, loading)
    fprintf('Creating 2D plot.\n');
    % Run advanced MVDR beamforming
    response_db = mvdr_beamforming(r, mic_positions, candidate_points, bin_freqs, c_0, loading);

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

    % Plot true source (red x) and estimated source (blue x)
    plot(source_x, source_y, 'rx', 'MarkerSize', 12, 'LineWidth', 2);
    plot(est_x, est_y, 'bx', 'MarkerSize', 12, 'LineWidth', 2);
    legend({'True Source','Estimated Source'}, 'Location','northeast');
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 14); % optional for report styling
end

% FUNCTION FIVE: 1D SCAN PLOTTING FOR MVDR
function plot_1dscan_mvdr(r, mic_positions, source_y, bin_freqs, c_0, loading)
    fprintf('Creating 1D plot for peak testing.\n');
    % Scan along y-axis at fixed x and z
    y_scan_line = linspace(-0.5, 0.5, 200);
    x_fixed = 0; 
    z_fixed = 0;

    % Build candidate points
    candidate_points_line = [x_fixed*ones(length(y_scan_line),1), ...
        y_scan_line(:), z_fixed*ones(length(y_scan_line),1)];

    % Run advanced MVDR beamforming along the line
    response_line = mvdr_beamforming(r, mic_positions, candidate_points_line, bin_freqs, c_0, loading);

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
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 14); % optional for report styling
end

