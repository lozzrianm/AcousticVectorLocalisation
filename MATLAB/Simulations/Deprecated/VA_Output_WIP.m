clc; clear all; close all;

%Generate array output using 3D sound field model
%For acoustic vector array of 1element each with four subcomponents

%Written by L Marshall 01/08/2025
%Edited by L Marshall 05/10/2025
%Edited by L Marshall 11/11/2025


%% DEFINE INPUT VARIABLES %%
% Signal Parameters
duration = 60; %Sample duration (s)
noise_level = 0.05; %Noise amplitude
num_sources = 1; %Number of signal sources

% More signal source parameters
% Source positions [x, y, z] (m)
source_positions = [
    0, -0.25, 0.0; %Source 1
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
   % 0.01; %Source 2 a
];

% Vector Sensor Array Parameters;
N_v = 1; %No. of vector sensors
d_y = 0.04; %Vector sensor y axis spacing (m)
delta = 0.04; %MEMS colocation spacing (m)

% Define geom centre of the array 
x_a = 2.0 + delta/2; % Geometric centre x-coordinate
z_a = 0; % Array centre z-coordinate
y_a = -d_y/2 + delta/2; % Geometric centre y-coordinate

% Physical Environment Parameters
c_o = 340; %Speed of sound (m/s)
rho_o = 1.02; %Air Density at STP

% Sampling Conditions
d_t = 0.0001; %Time step (s)
%T = ;
t_0 = 0.02; %Initial sampling time (s)
csv_filename = 'generatedsignal_avs.csv';

%% INITIALISE ARRARY GEOMETRY %%

% Calculate number of microphones
N_m = N_v * 4;  %4 MEMS mics per vector sensor

% Adjust array centre for geometry
if N_v == 1
    y_a_adjusted = y_a;
else
    % Centre the array
    y_a_adjusted = y_a - (N_v-1)*d_y/2;
end

%% DEFINE MIC POSITIONS %%
mic_positions = zeros(3, N_m);

for vs = 1:N_v
    idx = (vs-1)*4 + (1:4);
    
    % Centre of this vector sensor
    x_c = x_a;
    y_c = y_a_adjusted + (vs-1)*d_y;
    
    % Place four microphones in square configuration around centre
    mic_positions(:,idx(1)) = [x_c - delta/2; y_c - delta/2; z_a]; % p0: bottom-left
    mic_positions(:,idx(2)) = [x_c + delta/2; y_c - delta/2; z_a]; % p1: bottom-right
    mic_positions(:,idx(3)) = [x_c + delta/2; y_c + delta/2; z_a]; % p2: top-right
    mic_positions(:,idx(4)) = [x_c - delta/2; y_c + delta/2; z_a]; % p3: top-left
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

% Create header
header = {'Time'};
for vs = 1:N_v
    for mic = 1:4
        header{end+1} = sprintf('AVS%d_Mic%d_Monopole', vs, mic);
    end
end
for vs = 1:N_v
    for mic = 1:4
        header{end+1} = sprintf('AVS%d_Mic%d_Sinusoidal', vs, mic);
    end
end

% Save to CSV
results_table = array2table(results, 'VariableNames', header);
writetable(results_table, csv_filename);

fprintf('Saved signal data to: %s\n', csv_filename);

%% PLOT ARRAY GEOMETRY %%

figure('Name', 'Array Geometry');
hold on; grid on;

marker_size = 50;
marker_color = [0 0.447 0.741]; 

% Plot microphones grouped by vector sensor
for vs = 1:N_v
    mic_indices = (vs-1)*4 + (1:4);
    
    scatter(mic_positions(1, mic_indices), ...
            mic_positions(2, mic_indices), ...
            marker_size, marker_color, 'filled');
    
    % Label the first microphone of each sensor
    text(mic_positions(1, mic_indices(1)) - 0.05, ...
         mic_positions(2, mic_indices(1)) - 0.1, ...
         sprintf('VS %d', vs), 'FontSize', 10, 'FontName', 'Times New Roman');
end

% Plot all sources
source_colors = lines(num_sources);
for src = 1:num_sources
    scatter(source_positions(src, 1), ...
            source_positions(src, 2), ...
            80, source_colors(src,:), 'filled', 'MarkerEdgeColor', 'k');
    
    text(source_positions(src, 1) + 0.1, ...
         source_positions(src, 2), ...
         sprintf('S%d (%.0f Hz)', src, source_frequencies(src)), ...
         'FontSize', 10, 'FontName', 'Times New Roman');
end

xlabel('X (m)', 'FontName', 'Times New Roman', 'FontSize', 14);
ylabel('Y (m)', 'FontName', 'Times New Roman', 'FontSize', 14);

if num_sources == 1
    legend('Microphones', 'Source', 'Location', 'northwest');
else
    legend_entries = [{'Microphones'}, ...
        arrayfun(@(i) sprintf('Source %d', i), 1:num_sources, 'UniformOutput', false)];
    legend(legend_entries, 'Location', 'northwest');
end

set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
box on;
axis equal;


%% DISPLAY CONFIGURATION SUMMARY %%

fprintf('\n<strong>AVS Array Configuration</strong>\n');
fprintf('Number of AVS elements: %d\n', N_v);
fprintf('Microphones per AVS: 4 (for FD approximation)\n');
fprintf('Total microphones: %d\n', N_m);
fprintf('Finite difference spacing (delta): %.3f m\n', delta);
fprintf('AVS element spacing (d_y): %.3f m\n', d_y);

if num_sources == 1
    freq_display = source_frequencies(1);
    fprintf('Source frequency: %d Hz\n', freq_display);
    fprintf('Wavelength: %.3f m\n', c_o/freq_display);
    fprintf('Delta/wavelength ratio: %.3f\n', delta/(c_o/freq_display));
end

fprintf('\n<strong>Signal Sources: %d</strong>\n', num_sources);
for src = 1:num_sources
    fprintf('Source %d: (%.3f, %.3f, %.3f) m @ %.0f Hz\n', ...
        src, source_positions(src, :), source_frequencies(src));
end
