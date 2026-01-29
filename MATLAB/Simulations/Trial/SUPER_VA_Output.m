clc; close all; %clear all;

%Generate array output using 3D sound field model
%For acoustic vector array of 1element each with four subcomponents

%Written by L Marshall 27/11/2025

%% DEFINE INPUT VARIABLES %%

% Signal Parameters
duration = 30; %Sample duration (s)
noise_level = 0.05; %Noise amplitude
num_sources = 1; %Number of signal sources

% More signal source parameters
% Source positions [x, y, z] (m)
source_positions = [
    -0.7071, 0.7071, 0; %Source 1
    %4, -2, 0.0; %Source 2
];

% Source frequencies (Hz) - overridden in batch mode
if ~exist('batch_test_freq', 'var')
    source_frequencies = [
        1500; %Source 1 freq
        %2000 %Source 2 freq
    ];
else
    source_frequencies = batch_test_freq;
end

% Source amplitudes 
source_amplitudes = [
    0.01; %Source 1 a
   % 0.01; %Source 2 a
];

% Physical Environment Parameters
c_o = 340; %Speed of sound (m/s)
rho_o = 1.02; %Air Density at STP

% Vector Sensor Array(s) Parameters
N_a = 1; %Number of independent arrays
N_v = 1; %No. of vector sensors
d_y = 0.1; %Vector sensor y axis spacing (m)

% Delta - overridden in batch mode
if ~exist('batch_test_delta', 'var')
    delta = 0.042; %MEMS colocation spacing (m)
else
    delta = batch_test_delta;
end

% Array centre geom 1
x_a1 = 0;
z_a1 = 0; % Array centre z-coordinate
y_a1 = 0;
A1_coord = [x_a1; y_a1; z_a1];

% Array centre geom 2
lambda = c_o / source_frequencies(1, :);
a_spacing = sqrt(((5*lambda)^2)/2); %x and y spacing for 5*lambda distance
x_a2 = x_a1 - a_spacing; % Geometric centre x-coordinate
z_a2 = 0; % Array centre z-coordinate
y_a2 = y_a1 - a_spacing; % Geometric centre y-coordinate
A2_coord = [x_a2; y_a2; z_a2];

% Sampling Conditions
d_t = 0.0001; %Time step (s)

% CSV filename - overridden in batch mode
if ~exist('batch_csv_name', 'var')
    csv_filename = 'generatedsignal_avs.csv';
else
    csv_filename = batch_csv_name;
end

%% CHECK FOR BATCH MODE - SUPPRESS FIGURES %%

if exist('batch_test_freq', 'var')
    % Suppress figures in batch mode
    set(0, 'DefaultFigureVisible', 'off');
end
%% INITIALISE ARRARY GEOMETRY %%

% Calculate number of microphones
N_ma = N_v * 4; %4 MEMS per vector sensor in array
N_m = N_a * N_ma;  %total MEMS in system

% Adjust array centres for geometry
if N_v == 1
    y_a1_adjusted = y_a1;
    y_a2_adjusted = y_a2;
else
    % Centre each array along y-axis
    y_a1_adjusted = y_a1 - (N_v-1)*d_y/2;
    y_a2_adjusted = y_a2 - (N_v-1)*d_y/2;
end

%% DEFINE MIC POSITIONS %%
mic_positions = zeros(3, N_m);

% Store all array centres as columns
array_centres = [x_a1, x_a2;
                 y_a1_adjusted, y_a2_adjusted;
                 z_a1, z_a2];

% Loop through each independent array
for array = 1:N_a
    % Get this arrays centre coordinates
    centre = array_centres(:, array);
    
    % Loop through vector sensors
    for vs = 1:N_v
        % Calculate global microphone indices for this VS
        idx = (array-1)*N_ma + (vs-1)*4 + (1:4);
        
        % Centre of this vector sensor
        y_c = centre(2) + (vs-1)*d_y;
        
        % Place four microphones in square configuration around centre
        mic_positions(:, idx(1)) = [centre(1) - delta/2; y_c - delta/2; centre(3)]; %p0: bottom-left
        mic_positions(:, idx(2)) = [centre(1) + delta/2; y_c - delta/2; centre(3)]; %p1: bottom-right
        mic_positions(:, idx(3)) = [centre(1) + delta/2; y_c + delta/2; centre(3)]; %p2: top-right
        mic_positions(:, idx(4)) = [centre(1) - delta/2; y_c + delta/2; centre(3)]; %p3: top-left
    end
end

%% CREATE TIME VECTOR %%

num_samples = round(duration / d_t);
time = (0:num_samples - 1) * d_t;

fprintf('Generating %d samples over %.1f seconds\n', num_samples, duration);


%% INITIALISE SIGNAL ARRAYS %%

% Create storage for each source's contribution
% Dimensions: [time samples, microphones, sources]
p_monopole = zeros(length(time), N_m, num_sources);
p_sinusoidal = zeros(length(time), N_m, num_sources);


%% GENERATE SIGNALS FROM EACH SOURCE %%

fprintf('Generating signals from %d source(s)...\n', num_sources);

for src = 1:num_sources
    fprintf('  Processing source %d/%d (%.1f Hz)...', src, num_sources, ...
        source_frequencies(src));

    x_s = source_positions(src, 1);
    y_s = source_positions(src, 2);
    z_s = source_positions(src, 3);
    f_src = source_frequencies(src);
    Q_o = source_amplitudes(src);
    
    omega = f_src * 2 * pi; %Angular frequency
    
    % Generate signal at each microphone from this source
    for k = 1:num_samples
        t1 = time(k);
        
        for n = 1:N_m
            x_mic = mic_positions(1, n);
            y_mic = mic_positions(2, n);
            z_mic = mic_positions(3, n);
            
            % Calculate distance from source to microphone
            R_a = sqrt((x_mic - x_s)^2 + (y_mic - y_s)^2 + (z_mic - z_s)^2);
            
            % Account for propagation delay
            t = t1 - R_a/c_o;
            
            % Monopole derivative signal
            Q1_dot = Q_o * omega * cos(omega * t);
            p_monopole(k, n, src) = (rho_o / (4*pi*R_a)) * Q1_dot;
            
            % Sinusoidal signal
            Q2_dot = Q_o * sin(omega * t);
            p_sinusoidal(k, n, src) = (rho_o / (4*pi*R_a)) * Q2_dot;
        end
    end
    
    fprintf(' Done\n');
end


%% COMBINE SIGNALS FROM ALL SOURCES FOR DATA EXPORTING %%

% Sum contributions from all sources
p1 = sum(p_monopole, 3);     % Sum across source dimension
p2 = sum(p_sinusoidal, 3);

fprintf('Combined signals from all sources\n');


%% ADD SIMULATED NOISE %%

fprintf('Adding noise (level = %.3f)...\n', noise_level);

p1_noisy = p1 + noise_level * randn(size(p1));
p2_noisy = p2 + noise_level * randn(size(p2));


%% NORMALISE AND SAVE TO CSV %%

fprintf('Normalising and saving to CSV...\n');

% Normalise signals
generated_monopoles = p1_noisy ./ max(abs(p1_noisy), [], 1);
generated_sinusoidals = p2_noisy ./ max(abs(p2_noisy), [], 1);

% Prepare results matrix
results = [time(:), generated_monopoles, generated_sinusoidals];

% Create header for two independent arrays
header = {'Time'};

% Generate headers for monopole signals first, then sinusoidal
signal_types = {'Monopole', 'Sinusoidal'};

for sig_type = 1:2
    for array = 1:N_a
        for vs = 1:N_v
            for mic = 1:4
                header{end+1} = sprintf('Array%d_AVS%d_Mic%d_%s', ...
                    array, vs, mic, signal_types{sig_type});
            end
        end
    end
end

% Save to CSV
results_table = array2table(results, 'VariableNames', header);
writetable(results_table, csv_filename);

fprintf('Saved signal data to: %s\n', csv_filename);


%% PLOT ARRAY GEOMETRY %%

figure('Name', 'Array Geometry - Two Independent Arrays');
hold on; grid on;

marker_size = 50;

% Define distinct colours for each array
array_colours = [0 0.447 0.741; 0.850 0.325 0.098]; %Blue, Orange

% Plot microphones for each array
for array = 1:N_a
    colour = array_colours(array, :);
    
    for vs = 1:N_v
        % Calculate microphone indices for this VS
        mic_indices = (array-1)*N_ma + (vs-1)*4 + (1:4);
        
        scatter(mic_positions(1, mic_indices), ...
                mic_positions(2, mic_indices), ...
                marker_size, colour, 'filled');
        
        % Label the first microphone of the first sensor in each array
        if vs == 1
            text(mic_positions(1, mic_indices(1)) - 0.05, ...
                 mic_positions(2, mic_indices(1)) - 0.1, ...
                 sprintf('Array %d (VS %d)', array, vs), 'FontSize', 10, ...
                 'FontName', 'Times New Roman');
        else
            text(mic_positions(1, mic_indices(1)) - 0.05, ...
                 mic_positions(2, mic_indices(1)) - 0.1, ...
                 sprintf('VS %d', vs), 'FontSize', 10, ...
                 'FontName', 'Times New Roman');
        end
    end
end

% Plot all sources
source_colours = lines(num_sources);
for src = 1:num_sources
    scatter(source_positions(src, 1), ...
            source_positions(src, 2), ...
            80, source_colours(src,:), 'filled', 'MarkerEdgeColor', 'k');
    
    text(source_positions(src, 1) + 0.1, ...
         source_positions(src, 2), ...
         sprintf('S%d (%.0f Hz)', src, source_frequencies(src)), ...
         'FontSize', 10, 'FontName', 'Times New Roman');
end

xlabel('X (m)', 'FontName', 'Times New Roman', 'FontSize', 14);
ylabel('Y (m)', 'FontName', 'Times New Roman', 'FontSize', 14);

% Create legend entries
legend_entries = arrayfun(@(i) sprintf('Array %d', i), 1:N_a, 'UniformOutput', false);
if num_sources == 1
    legend_entries{end+1} = 'Source';
else
    source_entries = arrayfun(@(i) sprintf('Source %d', i), 1:num_sources, 'UniformOutput', false);
    legend_entries = [legend_entries, source_entries];
end
legend(legend_entries, 'Location', 'northwest');

set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
box on;

% Set appropriate axis limits with padding
all_x = [mic_positions(1,:), source_positions(:,1)'];
all_y = [mic_positions(2,:), source_positions(:,2)'];
x_range = max(all_x) - min(all_x);
y_range = max(all_y) - min(all_y);
padding = 0.15; %Padding as fraction of range (30%)

xlim([min(all_x) - padding*max(x_range,0.5), max(all_x) + padding*max(x_range,0.5)]);
ylim([min(all_y) - padding*max(y_range,0.5), max(all_y) + padding*max(y_range,0.5)]);

axis equal;


%% RESTORE FIGURE VISIBILITY %%

if exist('batch_test_freq', 'var')
    set(0, 'DefaultFigureVisible', 'on');
end

%% DISPLAY CONFIGURATION SUMMARY %%

fprintf('\n<strong>Dual Independent Array Configuration</strong>\n');
fprintf('Number of independent arrays: %d\n', N_a);
fprintf('Vector sensors per array: %d\n', N_v);
fprintf('Microphones per vector sensor: 4 (for FD approximation)\n');
fprintf('Microphones per array: %d\n', N_ma);
fprintf('Total microphones in system: %d\n', N_m);
fprintf('Finite difference spacing (delta): %.3f m\n', delta);

if N_v > 1
    fprintf('Vector sensor spacing (d_y): %.3f m\n', d_y);
end

% Array positioning information
array_separation = norm(A2_coord - A1_coord);
fprintf('\nArray 1 centre: (%.3f, %.3f, %.3f) m\n', x_a1, y_a1, z_a1);

if N_a > 1
    fprintf('Array 2 centre: (%.3f, %.3f, %.3f) m\n', x_a2, y_a2, z_a2);
    fprintf('Array separation distance: %.3f m\n', array_separation);
end

if num_sources == 1
    freq_display = source_frequencies(1);
    fprintf('\nSource frequency: %d Hz\n', freq_display);
    fprintf('Wavelength: %.3f m\n', c_o/freq_display);
    fprintf('Delta/wavelength ratio: %.3f\n', delta/(c_o/freq_display));
    fprintf('Array separation/wavelength ratio: %.3f\n', array_separation/(c_o/freq_display));
    
    if N_v > 1
        fprintf('VS spacing/wavelength ratio: %.3f\n', d_y/(c_o/freq_display));
    end
end

fprintf('\n<strong>Signal Sources: %d</strong>\n', num_sources);
for src = 1:num_sources
    fprintf('Source %d: (%.3f, %.3f, %.3f) m @ %.0f Hz\n', ...
        src, source_positions(src, :), source_frequencies(src));
end