clc; clear all; close all;

% L Marshall 15/09 
% Single Source Localisation via a Linear Microphone Array

% MVDR method of beamforming %
% Using broadband analysis with spectrum blocks %
% with eigenvalue decomposition, adaptive regularisation, and weighted MVDR

%SCRIPT ADAPTED ON 14/11 FOR EXPERIMENTAL CONFIGURATIONS%

%% DEFINE INPUT VARIABLES %%

% Input file parameters
wav_filename = '20251120-075620(UTC)-ULA4Sensor2011atorigin-0005916277.wav'; 
csv_filename = 'ULA4Sensor2011-arorigin.csv'; % Output CSV file name
convert_wav_to_csv = true; 

% Array Parameters
N_mics = 4; %number of receiving elements
d_y = 0.08; %sensor spacing (m)
x_a = 2.0; %array centre x coordinate (m)
z_a = 0; %array centre z coordinate (m)
y_a = 0; %array centre y coordinate (m)

% Physical Environment Parameters
c_0 = 340; %speed of sound (m/s)

% Sampling Characteristics
freq_limit = 3000; %max freq for analysis (Hz)
overlap = 0.5; %percentage overlap for sample window (%)
loading = 1e-4; %set regularisation parameter (loading var in 2020 script)

% Grid Search Parameters
x_scan_points = 200;
y_scan_points = 200;
x_margin = 3; %(m)
y_margin = 1; %(m)

% Visualisation Parameters
% Source position [x, y] (m) - for reference
source_x = 0;
source_y = 0;


%% CONVERT WAV TO CSV (if needed) %%

if convert_wav_to_csv
    fprintf('<strong>Converting WAV to CSV</strong>\n');
    
    % Check if WAV file exists
    if ~exist(wav_filename, 'file')
        error('WAV file not found: %s', wav_filename);
    end
    
    % Read WAV file
    fprintf('Reading WAV file: %s\n', wav_filename);
    [audio_data, sample_rate] = audioread(wav_filename);
    
    % Check number of channels
    [num_samples, num_channels] = size(audio_data);
    fprintf('  Samples: %d\n', num_samples);
    fprintf('  Channels: %d\n', num_channels);
    fprintf('  Sample rate: %d Hz\n', sample_rate);
    fprintf('  Duration: %.3f seconds\n', num_samples/sample_rate);
    
    % Create time vector
    time_vec = (0:num_samples-1)' / sample_rate;
    
    % Combine time and audio data
    csv_data = [time_vec, audio_data];
    
    % Write to CSV
    fprintf('Writing to CSV: %s\n', csv_filename);
    writematrix(csv_data, csv_filename);
    fprintf('Conversion complete!\n\n');
end

%% MICROPHONE CALIBRATION DATA PROCESSING %%

% Calibration Parameters
apply_calibration = true; %set to false to disable all corrections

% Signal Conditioning Parameters
apply_dc_removal = true; %remove DC offset from time-domain signals

% Measured absolute sensitivities (mV/Pa)
% Replace with your actual calibration data
absolute_sensitivity_mVPa = [
    49.0;  %mic 1
    52.0;  %mic 2
    50.0;  %mic 3
    47.0;  %mic 4
    % 50.7;  %mic 5
    % 50.0;  %mic 6
    % 48.7;  %mic 7
    % 51.6;  %mic 8
    % 49.5;  %mic 9
    % 50.4;  %mic 10
    % 50.1   %mic 11
];

% Normalise to the mean and calc relative sensitivities
mean_sensitivity = mean(absolute_sensitivity_mVPa);
mic_sensitivity = absolute_sensitivity_mVPa / mean_sensitivity; %column vector (11x1)

fprintf('<strong>Microphone Calibration Processing</strong>\n');
fprintf('Mean absolute sensitivity: %.2f mV/Pa\n', mean_sensitivity);
fprintf('\nConverting from Absolute to Relative Sensitivity:\n');
for mic = 1:length(mic_sensitivity)
    fprintf('  Mic %d: %.2f mV/Pa → %.4f (%.1f%% of mean)\n', ...
        mic, absolute_sensitivity_mVPa(mic), mic_sensitivity(mic), ...
        mic_sensitivity(mic) * 100);
end

% Phase offset corrections (radians)
% (+) vals mean microphone signal is ahead (leading)
% (-) vals mean microphone signal is behind (lagging)
% IMPORTANT: Must be column vector to match mic_sensitivity dimensions
phase_offset = [
    0;      %mic 1
    0;   %mic 2
   0;   %mic 3
    0;   %mic 4
   % -0.04;   %mic 5
   %  0.01;   %mic 6
   %  0.03;   %mic 7
   % -0.02;   %mic 8
   %  0.04;   %mic 9
   % -0.01;   %mic 10
   %  0       %mic 11
]; %column vector (11x1)

fprintf('\nPhase offsets: %.4f to %.4f rad (%.2f to %.2f deg)\n', ...
    min(phase_offset), max(phase_offset), ...
    rad2deg(min(phase_offset)), rad2deg(max(phase_offset)));
fprintf('Phase offset mean: %.4f rad (%.2f deg)\n\n', ...
    mean(phase_offset), rad2deg(mean(phase_offset)));

%% LOAD IN GENERATED SIGNAL FROM CSV %%

fprintf('Loading signal data from: %s\n', csv_filename);

if ~exist(csv_filename, 'file')
    error('Signal file not found: %s\nRun signal generation first.', csv_filename);
end

data = readmatrix(csv_filename);
time = data(:, 1); %time vector
[numRows, numCols] = size(data);

N = numRows; %number of samples
Nr = numCols-1; %number of receiving elements

fprintf('Loaded %d samples from %d microphones\n', N, Nr);

%% EXTRACT MICROPHONE SIGNALS %%

tx = data(:, 2:end).'; % transpose: rows = sensors, columns = time samples

fprintf('Signal matrix size: %d microphones x %d time samples\n', size(tx, 1), size(tx, 2));

%% SIGNAL CONDITIONING %%

if apply_dc_removal
    fprintf('Removing DC offset from signals...\n');
    tx = tx - mean(tx, 2); %subtract mean from each microphone signal
end

%% ARRAY GEOMETRY SET UP %%

n = linspace(-(Nr - 1) / 2, (Nr - 1) / 2, Nr); %number of sensors
mic_positions = [x_a * ones(1, Nr); d_y * n; zeros(1, Nr)]; %sensor positions
 
%% CONVERT TO FREQ DOMAIN %%

fprintf('Converting to frequency domain...\n');

d_t = time(2) - time(1); %time step (s)
F_s = 1 / d_t; %sampling frequency (Hz)
f = F_s * (0:floor(N/2)) / N; %frequency vector

% Limit frequency range
idx_limit = f <= freq_limit;

% Preallocate FFT matrix
tx_freq = zeros(Nr, sum(idx_limit));

% Plot frequency spectrum
figure('Name', 'Frequency Domain');
hold on;
for sensor = 1:Nr
    % FFT of each sensor signal
    Y = fft(tx(sensor, :));
    Y = Y(1:floor(N/2)+1); %keep positive freqs (complex)
    
    % Store only the useful freqs (<= freq_limit)
    tx_freq(sensor, :) = Y(idx_limit);
    
    % Plot magnitude spectrum 
    Y_mag = abs(Y) / N;
    Y_mag(2:end-1) = 2 * Y_mag(2:end-1);
    Y_dB = 20 * log10(Y_mag + eps);
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
z_scan = 0; %2D plane for now

[X_grid, Y_grid] = meshgrid(x_scan, y_scan);
Z_grid = zeros(size(X_grid));

candidate_points = [X_grid(:), Y_grid(:), Z_grid(:)]; %Nr x 3 matrix

% Grid resolution in x and y directions
x_res = x_scan(2) - x_scan(1);
y_res = y_scan(2) - y_scan(1);
grid_res = mean([x_res, y_res]);

fprintf('Search grid: %d x %d points, resolution: %.4f m\n', ...
    length(x_scan), length(y_scan), grid_res);

%% CREATING FREQ BINS FOR BROADBAND ANALYSIS %%

fprintf('Creating frequency bins...\n');

% Maximum bin width for req spatial resolution
max_binwidth = 1 / (8 * d_y * (Nr - 1) / c_0); %(Hz)
size_fft = floor(F_s / max_binwidth);

% Frequency spacing
delta_f = F_s / N; %frequency spacing per bin
fft_vec = F_s * (0:(size_fft - 1)) / size_fft; %freq vector

% Target frequencies
target_freqs = 0:max_binwidth:freq_limit;

% Map each target to nearest snapshot FFT bin
bin_index = zeros(size(target_freqs));
for k = 1:length(target_freqs)
    [~, bin_index(k)] = min(abs(fft_vec - target_freqs(k)));
end

% Generate bin indices within accepted range
bin_index = unique(bin_index); %remove duplicates
bin_index = bin_index(bin_index <= size_fft);
bin_freqs = fft_vec(bin_index); %actual bin freqs
num_bins = length(bin_index); %no. bins

fprintf('Using %d frequency bins (%.1f - %.1f Hz)\n', ...
    num_bins, min(bin_freqs), max(bin_freqs));

num_snaps = 4 * Nr; %no. snapshots

% Calc no. samples accounting for overlap
num_samps = (num_snaps + 1) * size_fft * overlap;

% Taper window
window = hanning(size_fft)';

%% CREATE SNAPSHOTS AND MAKE CSM %%

fprintf('Creating snapshots and cross-spectral matrix...\n');

snapshots = make_snapshots(tx, size_fft, overlap, window);
fprintf('  Created %d snapshots\n', size(snapshots, 3));

% Apply calibration corrections if enabled
if apply_calibration
    fprintf('  Applying calibration corrections...\n');
    snapshots = apply_calibration_correction(snapshots, mic_sensitivity, phase_offset);
end

% Max no. snapshots that signal fits
maxnum_snap = size(snapshots, 3); 

r = create_csm(snapshots, bin_index, Nr, num_bins, maxnum_snap);
fprintf('  Created %dx%dx%d cross-spectral matrix\n', size(r));

%% PERFORM BEAMFORMING %%

fprintf('\n<strong>Running MVDR Beamforming...</strong>\n');
fprintf('  Regularisation parameter: %.2e\n', loading);

% Run MVDR beamforming
response_db = mvdr_beamforming(r, mic_positions, candidate_points, bin_freqs, c_0, loading);

% 1D scan along y-axis at x=0 (for verification)
fprintf('\n1D Scan (x = 0):\n');
plot_1dscan_mvdr(r, mic_positions, source_y, bin_freqs, c_0, loading);

% 2D scan of full search area
fprintf('\n2D Scan:\n');
[est_x, est_y] = plot_2dscan_mvdr(r, mic_positions, candidate_points, bin_freqs, c_0, ...
    y_scan, x_scan, source_x, source_y, X_grid, Y_grid, loading);

%% PERFORMANCE METRICS %%

fprintf('\n<strong>Performance Metrics:</strong>\n');

est_pos = [est_x, est_y];
true_pos = [source_x, source_y];
errors = true_pos - est_pos;

squared_dist = sum(errors.^2);
MSE = squared_dist;
MSE_percent = (MSE / (grid_res^2)) * 100;
radial_error = sqrt(squared_dist);

% Angular error
y_sensor = 0; %array centre y coordinate
true_angle = atan2(source_y - y_sensor, source_x - x_a);
est_angle = atan2(est_y - y_sensor, est_x - x_a);
angular_error_deg = rad2deg(abs(true_angle - est_angle));
% Normalise to 180 deg
if angular_error_deg > 180
    angular_error_deg = mod(angular_error_deg + 180, 360) - 180;
end

fprintf('Estimated position: (%.3f, %.3f) m\n', est_x, est_y);
fprintf('Position error: dx = %.3f m, dy = %.3f m\n', errors(1), errors(2));
fprintf('Radial error: %.3f m\n', radial_error);
fprintf('Mean Squared Error: %.4f m²\n', MSE);
fprintf('MSE (normalised): %.2f%% of grid resolution\n', MSE_percent);
fprintf('Angular Error: %.4f deg\n', angular_error_deg);

%% FUNCTION DEFINITION %%

% FUNCTION: MAKE SNAPSHOTS OF GENERATED SIGNAL
% Make time domain snapshots of input signal
function snapshots = make_snapshots(tx, size_fft, overlap, window)
    % No. sample to move per snapshot
    step = round(overlap * size_fft);
    
    % No. of full snapshots that fit in signal with given step size
    num_snap = floor((size(tx, 2) - size_fft) / step) + 1; %no. snapshots 
    % Start indices for each snapshot
    start_idx = 1 + (0:(num_snap - 1)) * step;
    
    % Create index matrix (size_fft x num_snap)
    idx_matrix = start_idx + (0:size_fft - 1)'; %each column = indices for a snapshot
    
    % Extract snapshots (force reshape into 3D)
    snapshots = tx(:, idx_matrix);  
    snapshots = reshape(snapshots, size(tx, 1), size_fft, num_snap);
    
    % Apply window along snapshots dimension
    snapshots = snapshots .* reshape(window, [1, size_fft, 1]);
end


% FUNCTION: APPLY CALIBRATION CORRECTIONS
% Apply sensitivity and phase corrections in frequency domain
function snapshots_corrected = apply_calibration_correction(snapshots, mic_sensitivity, phase_offset)  
    N_mics = size(snapshots, 1);
    
    % FFT each snapshot to frequency domain
    snapshots_freq = fft(snapshots, [], 2); %(N_mics x size_fft x num_snap)
    
    % Create correction factors for each microphone
    % correction = (1/sensitivity) * exp(-i*phase_offset)
    % 1/sensitivity normalises amplitude
    % exp(-i*phase_offset) removes the phase error
    correction_factors = (1 ./ mic_sensitivity) .* exp(-1i * phase_offset);

    % Reshape to apply across all snapshots
    correction_factors = reshape(correction_factors, [N_mics, 1, 1]); 
    
    % Multiply each frequency bin by correction factor
    snapshots_freq_corrected = snapshots_freq .* correction_factors; 
    
    snapshots_corrected = ifft(snapshots_freq_corrected, [], 2);
    
    % Ensure real output
    snapshots_corrected = real(snapshots_corrected);
end


% FUNCTION: FFT AND CREATE CROSS-SPECTRAL MATRIX
% Convert each block to freq domain and create cross-spectral matrix
function r = create_csm(snapshots, bin_index, Nr, num_bins, maxnum_snap)
    % FFT along the time dimension for all snapshots
    fx = fft(snapshots, [], 2); %(Nr x size_fft x num_snap)    
    
    % Keep only the frequencies of interest
    fx = fx(:, bin_index, :); %(Nr x num_bins x num_snap)
    
    % Preallocate CSM
    r = zeros(Nr, Nr, num_bins);
    
    % Vectorised additions to CSM per snapshot
    for jf = 1:num_bins
        % Extract matrix for this frequency bin
        fx1 = squeeze(fx(:, jf, :)); %(Nr x num_snap)
    
        % Compute CSM: sum over snapshots of outer products
        r(:, :, jf) = fx1 * fx1'; %(Nr x Nr)
    end
    
    % Average CSM over snapshots of sampled signal
    r = r / maxnum_snap;
end


% FUNCTION: MVDR BEAMFORMING
function response_db = mvdr_beamforming(r, mic_positions, candidate_points, bin_freqs, c_0, loading)

    num_points = size(candidate_points, 1); %no. points to search
    num_bins = size(r, 3);

    % Preallocate arrays
    mvdr_responses = zeros(num_points, num_bins);
    cbf_responses = zeros(num_points, num_bins);

    for jf = 1:num_bins
        k = 2 * pi * bin_freqs(jf) / c_0; %wave no. for this freq bin
        rf = squeeze(r(:, :, jf)); %covariance matrix for this freq bin
        
        % Eigenvalue decomposition (improve stability)
        [ur, er] = eig(rf);
        erv = real(diag(er));
        
        % Safety check for negative eigenvalues 
        erv = max(erv, 1e-12);
        
        for n = 1:num_points
            l_pt = candidate_points(n, :).';
            
            % Calculate steering vector
            r_lm = sqrt(sum((mic_positions - l_pt).^2, 1));
            v1 = exp(-1i * k * r_lm).'; %steering vector 
            
            % Conventional DAS beamformer output 
            cbf_out = real(v1' * rf * v1);
            cbf_responses(n, jf) = cbf_out;
            
            % Adaptive regularisation based on CBF output
            lambda = loading * cbf_out; %lambda = diagonal loading factor
            if lambda == 0
                lambda = loading; %default to fixed loading
            end
            
            % Weighted MVDR using eigendecomposition
            rxv = (ur * diag(1 ./ (erv + lambda)) * ur') * v1;
            
            % Normalise to get weight vector
            denominator = v1' * rxv;
            if abs(denominator) < 1e-12
                mvdr_responses(n, jf) = 0;
            else
                vmvdr = rxv / denominator; %MVDR weight vector
                
                % Calc weighted output power
                mvdr_responses(n, jf) = abs(vmvdr' * rf * vmvdr);
            end
        end
    end

    % Broadband response by summing across frequencies
    responses_sum = sum(mvdr_responses, 2);

    % Convert to dB and normalise
    response_db = 10 * log10(responses_sum + eps);
    response_db = response_db - max(response_db);
end


% FUNCTION: 2D SCAN PLOTTING
function [est_x, est_y] = plot_2dscan_mvdr(r, mic_positions, candidate_points, ...
    bin_freqs, c_0, y_scan, x_scan, source_x, source_y, X_grid, Y_grid, loading)
    
    % Run MVDR beamforming
    response_db = mvdr_beamforming(r, mic_positions, candidate_points, bin_freqs, c_0, loading);

    % Reshape responses to a grid for plotting
    grid_response = reshape(response_db, length(y_scan), length(x_scan));

    % Plot spatial response map
    figure('Name', '2D MVDR Beamforming');
    imagesc(x_scan, y_scan, grid_response);
    axis xy;
    xlabel('X (m)', 'FontName', 'Times New Roman', 'FontSize', 14);
    ylabel('Y (m)', 'FontName', 'Times New Roman', 'FontSize', 14);
    colorbar;
    colormap('jet');
    hold on;

    % Find estimated source coordinates
    [max_val, max_idx] = max(grid_response(:));
    [y_idx, x_idx] = ind2sub(size(grid_response), max_idx);
    est_x = X_grid(y_idx, x_idx);
    est_y = Y_grid(y_idx, x_idx);

    fprintf('  Max response: %.2f dB at (%.3f, %.3f) m\n', max_val, est_x, est_y);

    % Plot true source (red x) and estimated source (blue x)
    plot(source_x, source_y, 'rx', 'MarkerSize', 12, 'LineWidth', 2, ...
        'DisplayName', 'True Source');
    plot(est_x, est_y, 'bx', 'MarkerSize', 12, 'LineWidth', 2, ...
        'DisplayName', 'Estimated Source');
    legend('Location', 'northeast');
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
end

% FUNCTION: 1D SCAN PLOTTING
function plot_1dscan_mvdr(r, mic_positions, source_y, bin_freqs, c_0, loading)
    
    % Scan along y-axis at fixed x and z
    y_scan_line = linspace(source_y - 0.5, source_y + 0.5, 200);
    x_fixed = 0; 
    z_fixed = 0;

    % Build candidate points
    candidate_points_line = [x_fixed * ones(length(y_scan_line), 1), ...
        y_scan_line(:), z_fixed * ones(length(y_scan_line), 1)];

    % Run MVDR beamforming along the line
    response_line = mvdr_beamforming(r, mic_positions, candidate_points_line, bin_freqs, c_0, loading);

    % Find maximum response
    [max_response, max_idx] = max(response_line);
    y_est = y_scan_line(max_idx);

    fprintf('  Max response: %.2f dB at y = %.3f m\n', max_response, y_est);
    
    % Plot the line response
    figure('Name', '1D MVDR Scan');
    plot(y_scan_line, response_line, 'LineWidth', 2);
    xlabel('y position (m)', 'FontName', 'Times New Roman', 'FontSize', 14);
    ylabel('Beamformer Response (dB)', 'FontName', 'Times New Roman', 'FontSize', 14);
    grid on;
    hold on;

    % Plot true source and estimated source
    plot(source_y, max(response_line), 'kx', 'MarkerSize', 10, 'LineWidth', 1.5); 
    xline(source_y, 'r--', 'LineWidth', 1.5, 'DisplayName', 'True Source');
    legend('Beamformer Response', 'True Source');
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
end
