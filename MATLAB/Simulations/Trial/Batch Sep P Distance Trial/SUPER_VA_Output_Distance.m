clc; close all; %clear all;

% L Marshall 27/11/2025
% Generate array output using 3D sound field model
% Modified for Distance Sweep - accepts batch_test_source_x and batch_test_source_y

%% DEFINE INPUT VARIABLES %%

% Signal Parameters
duration = 30;
noise_level = 0.05;
num_sources = 1;

% Source position - from batch script if available
if exist('batch_test_source_x', 'var') && exist('batch_test_source_y', 'var')
    source_positions = [batch_test_source_x, batch_test_source_y, 0];
else
    source_positions = [-0.7071, 0.7071, 0];
end

% Source frequencies (Hz)
if ~exist('batch_test_freq', 'var')
    source_frequencies = [1500];
else
    source_frequencies = batch_test_freq;
end

source_amplitudes = [0.01];

% Physical Environment Parameters
c_o = 340;
rho_o = 1.02;

% Vector Sensor Array Parameters
N_a = 1;
N_v = 1;
d_y = 0.1;

if ~exist('batch_test_delta', 'var')
    delta = 0.042;
else
    delta = batch_test_delta;
end

% Array centre
x_a1 = 0;
z_a1 = 0;
y_a1 = 0;
A1_coord = [x_a1; y_a1; z_a1];

lambda = c_o / source_frequencies(1, :);
a_spacing = sqrt(((5*lambda)^2)/2);
x_a2 = x_a1 - a_spacing;
z_a2 = 0;
y_a2 = y_a1 - a_spacing;
A2_coord = [x_a2; y_a2; z_a2];

% Sampling Conditions
d_t = 0.0001;

% CSV filename
if ~exist('batch_csv_name', 'var')
    csv_filename = 'generatedsignal_avs.csv';
else
    csv_filename = batch_csv_name;
end

%% CHECK FOR BATCH MODE - SUPPRESS FIGURES %%

if exist('batch_test_freq', 'var')
    set(0, 'DefaultFigureVisible', 'off');
end

%% INITIALISE ARRAY GEOMETRY %%

N_ma = N_v * 4;
N_m = N_a * N_ma;

if N_v == 1
    y_a1_adjusted = y_a1;
    y_a2_adjusted = y_a2;
else
    y_a1_adjusted = y_a1 - (N_v-1)*d_y/2;
    y_a2_adjusted = y_a2 - (N_v-1)*d_y/2;
end

%% DEFINE MIC POSITIONS %%

mic_positions = zeros(3, N_m);

array_centres = [x_a1, x_a2;
                 y_a1_adjusted, y_a2_adjusted;
                 z_a1, z_a2];

for array = 1:N_a
    centre = array_centres(:, array);
    
    for vs = 1:N_v
        idx = (array-1)*N_ma + (vs-1)*4 + (1:4);
        
        y_c = centre(2) + (vs-1)*d_y;
        
        mic_positions(:, idx(1)) = [centre(1) - delta/2; y_c - delta/2; centre(3)];
        mic_positions(:, idx(2)) = [centre(1) + delta/2; y_c - delta/2; centre(3)];
        mic_positions(:, idx(3)) = [centre(1) + delta/2; y_c + delta/2; centre(3)];
        mic_positions(:, idx(4)) = [centre(1) - delta/2; y_c + delta/2; centre(3)];
    end
end

%% CREATE TIME VECTOR %%

num_samples = round(duration / d_t);
time = (0:num_samples - 1) * d_t;

if ~exist('batch_test_freq', 'var')
    fprintf('Generating %d samples over %.1f seconds\n', num_samples, duration);
end

%% INITIALISE SIGNAL ARRAYS %%

p_monopole = zeros(length(time), N_m, num_sources);
p_sinusoidal = zeros(length(time), N_m, num_sources);

%% GENERATE SIGNALS FROM EACH SOURCE %%

if ~exist('batch_test_freq', 'var')
    fprintf('Generating signals from %d source(s)...\n', num_sources);
end

for src = 1:num_sources
    if ~exist('batch_test_freq', 'var')
        fprintf('  Processing source %d/%d (%.1f Hz)...', src, num_sources, ...
            source_frequencies(src));
    end

    x_s = source_positions(src, 1);
    y_s = source_positions(src, 2);
    z_s = source_positions(src, 3);
    f_src = source_frequencies(src);
    Q_o = source_amplitudes(src);
    
    omega = f_src * 2 * pi;
    
    for k = 1:num_samples
        t1 = time(k);
        
        for n = 1:N_m
            x_mic = mic_positions(1, n);
            y_mic = mic_positions(2, n);
            z_mic = mic_positions(3, n);
            
            R_a = sqrt((x_mic - x_s)^2 + (y_mic - y_s)^2 + (z_mic - z_s)^2);
            
            t = t1 - R_a/c_o;
            
            Q1_dot = Q_o * omega * cos(omega * t);
            p_monopole(k, n, src) = (rho_o / (4*pi*R_a)) * Q1_dot;
            
            Q2_dot = Q_o * sin(omega * t);
            p_sinusoidal(k, n, src) = (rho_o / (4*pi*R_a)) * Q2_dot;
        end
    end
    
    if ~exist('batch_test_freq', 'var')
        fprintf(' Done\n');
    end
end

%% COMBINE SIGNALS FROM ALL SOURCES %%

p1 = sum(p_monopole, 3);
p2 = sum(p_sinusoidal, 3);

if ~exist('batch_test_freq', 'var')
    fprintf('Combined signals from all sources\n');
end

%% ADD SIMULATED NOISE %%

if ~exist('batch_test_freq', 'var')
    fprintf('Adding noise (level = %.3f)...\n', noise_level);
end

p1_noisy = p1 + noise_level * randn(size(p1));
p2_noisy = p2 + noise_level * randn(size(p2));

%% NORMALISE AND SAVE TO CSV %%

if ~exist('batch_test_freq', 'var')
    fprintf('Normalising and saving to CSV...\n');
end

generated_monopoles = p1_noisy ./ max(abs(p1_noisy), [], 1);
generated_sinusoidals = p2_noisy ./ max(abs(p2_noisy), [], 1);

results = [time(:), generated_monopoles, generated_sinusoidals];

header = {'Time'};

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

results_table = array2table(results, 'VariableNames', header);
writetable(results_table, csv_filename);

if ~exist('batch_test_freq', 'var')
    fprintf('Saved signal data to: %s\n', csv_filename);
end

%% RESTORE FIGURE VISIBILITY %%

if exist('batch_test_freq', 'var')
    set(0, 'DefaultFigureVisible', 'on');
end