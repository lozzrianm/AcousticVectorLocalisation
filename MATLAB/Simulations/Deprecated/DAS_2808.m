clc; clear all; close all;

% L Marshall 28/08 
% Inverse Point Localisation of Linear Microphone Array

% Delay-and-sum method of beamforming %

%% Load in Generated Signal from CSV %%
% Data generated in 'Array_Output_2307.m'

data = readmatrix('generatedsignal.csv');
time = data(:,1); %Time vector
[numRows, numCols] = size(data);

N = numRows; %Number of samples
Nr = (numCols - 1)/2; %Number of receiving elements 

% Preallocate transmitted signal array
tx = zeros(N, Nr); 

% % Only capture monopole signals
% for sensor = 2:(Nr+1)
%     tx(:,sensor-1) = data(:, sensor);
% end

for mic = 1:Nr
    % Monopole Signals
    monopole_col = mic + 1; %Columns 2-12
    % Modulated Signals
    modulated_col = mic + 12; %Columns 13-23

    tx(:, mic) = data(:, monopole_col) + data(:, modulated_col);
end

tx = tx.'; %Transpose signals

%% Initialising Array and Signal Parameters %%

d_t = time(2) - time(1); %Time step (s)
sample_rate = 1/d_t; %Sampling rate

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

F_s = 1/(time(2)-time(1)); %Sampling freq (Hz)
f = F_s*(0:floor(N/2))/N; %Freq vector

% Limit frequency range to 0–1000 Hz
freq_limit = 1000;
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
title('Transmitted Signal in Frequency Domain');
xlim([0 freq_limit]);
legend show;
grid on;

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


%% Delay-and-Sum Beamforming %%
f_limited = f(idx_limit); %Frequency vector match tx_freq

% Find all viable peaks in spectrum
target_freqs = search_signalpeaks(f_limited, tx_freq);

% Loop through viable peaks for testing
for bin = 1:length(target_freqs)
  
  f0 = target_freqs(bin); %Set the peak freq
  fprintf('<strong>DAS Beamforming Response at %.1f Hz </strong>\n', f0)

  [~, f_peak_idx] = find(f_limited == f0); %Correlate to the index

  % New wave number
  k = (2*pi*f0)/c_0;
    
  % Select only the peak bin across all sensors
  tx_freq_peak = tx_freq(:, f_peak_idx);

  % Run DAS beamforming
  response_db = das_beamforming(tx_freq_peak, mic_positions, candidate_points, k);
  % Plot the 1D scan of beamformed response along x=0 to check
  plot_1dscan(tx_freq_peak, mic_positions, source_y, k, f0);
  % Plot the 2D scan of beamformed response
  plot_2dscan(response_db, y_scan, x_scan, f0, source_x, source_y, X_grid, Y_grid);

end


%% Function Definition %%

% FUNCTION 1: DAS_BEAMFORMING %
function response_db = das_beamforming(tx_freq, mic_positions, candidate_points, k)
num_points = size(candidate_points, 1);

responses = zeros(num_points, 1);

for n = 1:num_points
    l_pt = candidate_points(n,:).'; %Accesses vector [x y z]

    % Find distance between scanned point l and each sensor m
    r_lm = sqrt(sum((mic_positions-l_pt).^2,1));
    % Find distance between scanned point l and array centre
    %a_centre = mean(mic_positions,2);
    %r_l0 = norm(l_pt - a_centre);

    % Apply steering vector with phase shifts
    steer_vec = exp(-1i*k*r_lm).';

    beamformed_signal = abs(sum(conj(steer_vec) .* tx_freq))^2;

    responses(n) = beamformed_signal;

end

% Convert to dB scale
response_db = 20*log10(responses);
response_db = response_db - max(response_db);

end


% FUNCTION 2: FINDING INTERESTING PEAKS %
function target_freqs = search_signalpeaks(f_limited, tx_freq)

% Spectrum of the first sensor 
Y_mag = abs(tx_freq(1,:)); 
Y_dB  = 20*log10(Y_mag); %Convert to dB

% Put all points above 10 dB into a candidate list
candidate_idx = find(Y_dB > 10);
candidate_freqs = f_limited(candidate_idx);
candidate_mags  = Y_dB(candidate_idx);
all_peaks = [candidate_freqs(:), candidate_mags(:)];

% Initialise selected peaks array
selected_peaks = [];

% Loop through candidates and apply spacing check
for i = 1:size(all_peaks,1)
    f_i = all_peaks(i,1);
    mag_i = all_peaks(i,2);

    if isempty(selected_peaks)
        selected_peaks = all_peaks(i,:);
        continue
    end

    % Check frequency difference to existing stored peaks
    freq_diffs = abs(selected_peaks(:,1) - f_i);
    [min_diff, idx] = min(freq_diffs);

    if min_diff < 100  %too close
        % Keep the one with larger magnitude
        if mag_i > selected_peaks(idx,2)
            selected_peaks(idx,:) = all_peaks(i,:);
        end
    else
        % Add as a new peak
        selected_peaks = [selected_peaks; all_peaks(i,:)];
    end
end

% Display the selected peaks
fprintf('Detected peaks >10 dB and spaced >100 Hz:\n');
for i = 1:size(selected_peaks,1)
    fprintf('Peak at %.2f Hz, amplitude = %.2f dB\n', selected_peaks(i,1), selected_peaks(i,2));
end

% Use these frequencies for DAS beamforming
target_freqs = selected_peaks(:,1);

end



% FUNCTION 3: 2D SCAN PLOTTING %
function plot_2dscan(response_db, y_scan, x_scan, f0, source_x, source_y, X_grid, Y_grid)

    % Reshape responses to a grid for plotting
    grid_response = reshape(response_db, length(y_scan), length(x_scan));
    
    % Plot spatial response map
    figure;
    imagesc(x_scan, y_scan, grid_response);
    axis xy;
    xlabel('X (m)');
    ylabel('Y (m)');
    title(sprintf('DAS Beamforming Response at %.1f Hz', f0));
    colorbar;
    colormap('jet');
    
    hold on;
    
    % Find estimated source coords
    [max_val, max_idx] = max(grid_response(:)); %Global max
    [y_idx, x_idx] = ind2sub(size(grid_response), max_idx);
    est_x = X_grid(y_idx, x_idx);
    est_y = Y_grid(y_idx, x_idx);
    
    % For debugging
    fprintf('2D global max: %.2f dB at x = %.3f m, y = %.3f m\n', max_val, est_x, est_y);
    
    % Plot true source (red x)
    plot(source_x, source_y, 'rx', 'MarkerSize', 12, 'LineWidth', 2);
    hold on; 
    
    % Plot estimated source (blue x)
    plot(est_x, est_y, 'bx', 'MarkerSize', 12, 'LineWidth', 2);
    legend({'True Source','Estimated Source'}, 'Location','northeast');

end


% FUNCTION 4: 1D SCAN PLOTTING %
function plot_1dscan(tx_freq_peak, mic_positions, source_y, k, f0)

    y_scan_line = linspace(-0.5, 0.5, 200);
    x_fixed = 0; 
    z_fixed = 0;
    
    % Build candidate points (x=0 fixed, z=0 fixed, y varies)
    candidate_points_line = [x_fixed*ones(length(y_scan_line),1), ...
                             y_scan_line(:), ...
                             z_fixed*ones(length(y_scan_line),1)];
    
    % Run DAS beamforming for this line
    response_line = das_beamforming(tx_freq_peak, mic_positions, candidate_points_line, k);
    
    % Find maximum response
    [max_response, max_idx] = max(response_line);
    y_est = y_scan_line(max_idx);
    
    % Print results to command window
    fprintf('1D scan: max response = %.2f dB at y = %.3f m\n', max_response, y_est);
    
    % Plot the line response
    figure;
    plot(y_scan_line, response_line, 'LineWidth', 2);
    xlabel('y position (m)');
    ylabel('Beamformer Response (dB)');
    title(sprintf('1D Beamforming along y-axis at %.1f Hz', f0));
    grid on;
    
    % Label true source and estimated source
    hold on;
    plot(source_y, max(response_line), 'kx'); 
    xline(source_y, 'r--', 'LineWidth', 1.5, 'DisplayName','True Source');
    legend('Beamformer Response','True Source');

end



