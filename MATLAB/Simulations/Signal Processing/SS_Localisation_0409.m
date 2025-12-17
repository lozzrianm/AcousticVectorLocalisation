clc; clear all; close all;

% L Marshall 04/09 
% SS_Localisation_0409.m
% Single Source Localisation via a Linear Microphone Array

% Delay-and-sum method of beamforming %
% Using broadband analysis with spectrum blocks %

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
    sinusoidal_col = mic + (Nr+1); %Columns 13-23

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
source_positions = [source_x, source_y,0];

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
%title('Transmitted Signal in Frequency Domain');
xlim([0 freq_limit]);
legend show;
grid on;

%% Defining the Search Area %%
% Note: need to make it around the array center
x_margin = 3;
y_margin = 2;

x_scan = linspace(x_a-x_margin, x_a+x_margin, 200); 
y_scan = linspace(y_a-y_margin, y_a+y_margin, 200);

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


%% Perform Delay-and-Sum Beamforming %%

fprintf('<strong>DAS Beamforming Response in Broadband Spectrum</strong>\n')

% Run DAS beamforming
response_db = das_beamforming_csm(r, mic_positions, candidate_points, bin_freqs, c_0);
% Plot the 1D scan of beamformed response along x=0 to check
plot_1dscan(r, mic_positions, source_y, bin_freqs, c_0)
% Plot the 2D scan of beamformed response
[est_x, est_y] = plot_2dscan(r, mic_positions, candidate_points, bin_freqs, c_0, ...
    y_scan, x_scan, source_x, source_y, X_grid, Y_grid);

est_pos = [est_x est_y];
true_pos = [source_x source_y];
errors = true_pos - est_pos;
squared_dist = sum(errors.^2, 2); %Squared distance per source
MSE = mean(squared_dist); %Average across sources
% Normalised MSE w.r.t grid resolution
MSE_percent = (MSE / (grid_res^2)) * 100;

fprintf('Mean Squared Error %.4f m^2\n', MSE);
fprintf('Mean Squared Error (normalised): %.2f %% of ave. grid resolution\n', MSE_percent);

%% CALCULATE DETAILED PERFORMANCE METRICS %%

fprintf('\n<strong>PERFORMANCE METRICS:</strong>\n');

% Positions
true_x = source_x;
true_y = source_y;

fprintf('X Estimated position (m):         %.4f\n', est_x);
fprintf('Y Estimated Position (m):         %.4f\n', est_y);

% Position errors
dx = true_x - est_x;
dy = true_y - est_y;

fprintf('X Position Error, dx (m):         %.4f\n', dx);
fprintf('Y Position Error, dy (m):         %.4f\n', dy);

% Radial error
radial_error = sqrt(dx^2 + dy^2);
fprintf('Radial Error (m):                 %.4f\n', radial_error);

% Mean Squared Error
MSE = dx^2 + dy^2;  % Same as radial_error^2
fprintf('Mean Squared Error (m²):          %.4f\n', MSE);

% Normalised MSE
MSE_percent = (MSE / (grid_res^2)) * 100;
fprintf('Normalised as %% of Grid Resolu:  %.2f%%\n', MSE_percent);

% Angular error
% Reference point is array centre
array_centre_x = x_a;
array_centre_y = y_a;

true_angle = atan2(true_y - array_centre_y, true_x - array_centre_x);
est_angle = atan2(est_y - array_centre_y, est_x - array_centre_x);
angular_error_deg = rad2deg(abs(true_angle - est_angle));

% Normalise to 180 degrees
if angular_error_deg > 180
    angular_error_deg = 360 - angular_error_deg;
end

fprintf('Absolute Angular Error (deg):     %.3f\n', angular_error_deg);

% Peak response at estimated position
grid_response = reshape(response_db, length(y_scan), length(x_scan));
[~, est_idx] = min(abs(X_grid(:) - est_x) + abs(Y_grid(:) - est_y));
peak_response = grid_response(est_idx);

fprintf('Peak Response (dB):               %.2f\n', peak_response);
fprintf('\n');


%% BEAM PATTERN ANALYSIS %%
fprintf('\n<strong>GENERATING BEAM PATTERNS</strong>\n');

array_centre = [x_a; y_a];
plot_beam_pattern(response_db, X_grid, Y_grid, ...
    array_centre, "Array Beam Pattern", source_positions(:,1:2));

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
% Based on 2020 code by Dr.Bao
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

% FUNCTION THREE: DAS BEAMFORMING %
function response_db = das_beamforming_csm(r, mic_positions, candidate_points, bin_freqs, c_0)
% Notes on sizing:
% r: Nr x Nr x num_bins cross-spectral matrix
% mic_positions: 3 x Nr
% candidate_points: N_points x 3

num_points = size(candidate_points,1); %No. points to search
num_bins = size(r,3);

responses = zeros(num_points, num_bins);

for jf = 1:num_bins
    k = 2*pi*bin_freqs(jf)/c_0; %Wave no. for this freq bin

    for n = 1:num_points
        l_pt = candidate_points(n,:).'; %Accesses vector [x y z]
        % Find distance between scanned point l and each sensor m
        r_lm = sqrt(sum((mic_positions-l_pt).^2,1));

        % Find distance between scanned point l and array centre
        %a_centre = mean(mic_positions,2);
        %r_l0 = norm(l_pt - a_centre);

        % Apply steering vector with phase shifts
        steer_vec = exp(-1i*k*r_lm).';

        responses(n,jf) = abs(steer_vec'*r(:,:,jf)*steer_vec);
        
    end
end

% Sum over frequencies for broadband
responses_sum = sum(responses,2);

% Convert to dB and normalise
response_db = 20*log10(responses_sum);
response_db = response_db-max(response_db);

end


% FUNCTION FOUR: 2D SCAN PLOTTING %
function [est_x, est_y] = plot_2dscan(r, mic_positions, candidate_points, ...
    bin_freqs, c_0, y_scan, x_scan, source_x, source_y, X_grid, Y_grid)

    % Run DAS beamforming using CSM
    response_db = das_beamforming_csm(r, mic_positions, candidate_points, bin_freqs, c_0);

    % Reshape responses to a grid for plotting
    grid_response = reshape(response_db, length(y_scan), length(x_scan));

    % Plot spatial response map
    figure;
    imagesc(x_scan, y_scan, grid_response);
    axis xy;
    xlabel('X (m)');
    ylabel('Y (m)');
    %title('DAS Beamforming Response (Broadband)');
    colorbar;
    colormap('jet');
    hold on;

    % % Plot array sensor positions
    % plot(mic_positions(1,:), mic_positions(2,:), 'ko', 'MarkerFaceColor','k', ...
    %     'MarkerSize', 6);
    % hold on;

    % Find estimated source coords
    [max_val, max_idx] = max(grid_response(:)); %Global max
    [y_idx, x_idx] = ind2sub(size(grid_response), max_idx);
    est_x = X_grid(y_idx, x_idx);
    est_y = Y_grid(y_idx, x_idx);

    % For debugging
    fprintf('2D global max: %.2f dB at x = %.3f m, y = %.3f m\n', max_val, est_x, est_y);

    % Plot true source (red x) and estimated source (blue x)
    plot(source_x, source_y, 'rx', 'MarkerSize', 12, 'LineWidth', 2);
    plot(est_x, est_y, 'bx', 'MarkerSize', 12, 'LineWidth', 2);
    % legend({'Sensors','True Source','Estimated Source'}, 'Location','northeast');
    legend({'True Source','Estimated Source'}, 'Location','northeast');
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 14); % optional for report styling
end


% FUNCTION FIVE: 1D SCAN PLOTTING %
function plot_1dscan(r, mic_positions, source_y, bin_freqs, c_0)
    % Scan along y-axis at fixed x and z
    y_scan_line = linspace(-0.5, 0.5, 200);
    x_fixed = 0; 
    z_fixed = 0;

    % Build candidate points (x=0 fixed, z=0 fixed, y varies)
    candidate_points_line = [x_fixed*ones(length(y_scan_line),1), ...
        y_scan_line(:), z_fixed*ones(length(y_scan_line),1)];

    % Run DAS beamforming along the line
    response_line = das_beamforming_csm(r, mic_positions, candidate_points_line, bin_freqs, c_0);

    % Find maximum response
    [max_response, max_idx] = max(response_line);
    y_est = y_scan_line(max_idx);

    % For debugging
    fprintf('1D scan  at x=0: max response = %.2f dB at y = %.3f m\n', max_response, y_est);

    % Plot the line response
    figure;
    plot(y_scan_line, response_line, 'LineWidth', 2);
    xlabel('y position (m)');
    ylabel('Beamformer Response (dB)');
    %title('1D Beamforming along y-axis in broadband');
    grid on;
    hold on;

    % Plot true source and estimated source
    plot(source_y, max(response_line), 'kx', 'MarkerSize', 10, 'LineWidth', 1.5); 
    xline(source_y, 'r--', 'LineWidth', 1.5, 'DisplayName','True Source');
    legend('Beamformer Response','True Source');
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 14); % optional for report styling
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
            if source_angle_rad < 0
                source_angle_rad = source_angle_rad + 2*pi;
            end
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
