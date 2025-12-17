clc; clear all; close all;

% Generate array output using 3D sound field model
% Adapted from code by J Pan 24/6/2025
% Written by L Marshall 01/08/2025
% Edited by L Marshall 04/09/2025
% Edited for organisation by L Marshall 14/11/2025

%% DEFINE INPUT VARIABLES %%

% Signal Parameters
duration = 60; %sample duration (s)
noise_level = 0.05; %noise amplitude
Q_o = 0.01; %source amplitude
f = 1000; %source frequency (Hz)

% Source Position
x_s = 0; %source x coordinate (m)
y_s = -0.25; %source y coordinate (m)
z_s = 0; %source z coordinate (m)

% Array Parameters
N_mics = 11; %number of microphones
d_y = 0.08; %spacing between array elements (m)
x_a = 2; %array centre x coordinate (m)
z_a = 0; %array centre z coordinate (m)

% Physical Environment Parameters
c_o = 340; %speed of sound (m/s)
roh_o = 1.02; %air density at STP

% Sampling Conditions
d_t = 0.0001; %time step (s)
csv_filename = 'generatedsignal.csv';

%% INITIALISE SIGNAL CHARACTERISTICS %%

omega = f * 2 * pi; %angular frequency

%% CREATE TIME VECTOR %%

num_samples = round(duration / d_t);
time = (0:num_samples - 1) * d_t;

fprintf('Generating %d samples over %.1f seconds\n', num_samples, duration);

%% INITIALISE PRESSURE ARRAYS %%

p1 = zeros(length(time), N_mics); %monopole derivative
p2 = zeros(length(time), N_mics); %sinusoidal

%% GENERATE SIGNALS %%

fprintf('Generating signals...\n');

for k = 1:num_samples
    t1 = time(k);
    
    % Array output
    for n = 1:N_mics
        y_a = (n - 6) * d_y;
        R_a = sqrt((x_a - x_s)^2 + (y_a - y_s)^2 + (z_a - z_s)^2);
        t = t1 - R_a / c_o;
        
        % Sinusoidal source derivative signals
        Q1_dot = Q_o * omega * cos(omega * t); %monopole derivative
        Q2_dot = Q_o * sin(omega * t); %sinusoidal source
        
        p1(k, n) = (roh_o / (4 * pi * R_a)) * Q1_dot;
        p2(k, n) = (roh_o / (4 * pi * R_a)) * Q2_dot;
    end
end

fprintf('Signal generation complete\n');

%% ADD SIMULATED NOISE %%

fprintf('Adding noise (level = %.3f)...\n', noise_level);

p1_noisy = p1 + noise_level * randn(size(p1));
p2_noisy = p2 + noise_level * randn(size(p2));

%% PLOT TIME DOMAIN SIGNALS %%

figure('Name', 'Array Output Time Domain');
plot(time, p1_noisy(:,1) / max(p1_noisy(:,1)), 'k-', ...
     time, p2_noisy(:,1) / max(p2_noisy(:,1)), 'r-', 'LineWidth', 1);
ylabel('Array Output (V)', 'FontName', 'Times New Roman', 'FontSize', 14);
xlabel('Time (Seconds)', 'FontName', 'Times New Roman', 'FontSize', 14);
legend('Monopole Derivative', 'Sinusoidal Source');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
xlim([0 0.01]);
ylim([-1.5 1.5]);

%% NORMALISE AND SAVE TO CSV %%

fprintf('Normalising and saving to CSV...\n');

generated_monopoles = p1_noisy ./ max(p1_noisy, [], 1);
generated_sinusoidals = p2_noisy ./ max(p2_noisy, [], 1);

results = [time(:), generated_monopoles, generated_sinusoidals];

% Save data to file
writematrix(results, csv_filename);

fprintf('Saved signal data to: %s\n', csv_filename);

%% CALCULATE DIRECTION OF ARRIVAL %%

% Array centre position
x_array = x_a;
y_array = 0;

% Vector from source to array centre
dx = x_array - x_s;
dy = y_array - y_s;

% Direction of arrival (DoA) angle in radians and degrees
theta_rad = atan2(-dy, dx); %negative dy because y-axis increases upward
theta_deg = rad2deg(theta_rad);

fprintf('Direction of Arrival (DoA): %.2f degrees\n', theta_deg);

%% PLOT ARRAY GEOMETRY %%

figure('Name', 'Array Geometry');
hold on;
grid on;
axis equal;

% Microphone positions
y_positions = (1:N_mics) - 2.5;
y_positions = y_positions * d_y;
x_positions = x_a * ones(1, N_mics);

% Plot microphones
scatter(x_positions, y_positions, 80, 'b', 'filled');

% Plot source
scatter(x_s, y_s, 100, 'r', 'filled');

% Plot lines from source to each microphone
for n = 1:N_mics
    plot([x_s, x_positions(n)], [y_s, y_positions(n)], ...
        '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5);
end

% Axes labels and formatting
xlabel('X (m)', 'FontName', 'Times New Roman', 'FontSize', 14);
ylabel('Y (m)', 'FontName', 'Times New Roman', 'FontSize', 14);
legend('Microphones', 'Source', 'Location', 'northwest');
xlim([-0.3, 2.1]);
ylim([-0.4, 0.4]);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
box on;

%% DISPLAY CONFIGURATION SUMMARY %%

fprintf('\n<strong>Linear Array Configuration</strong>\n');
fprintf('Number of microphones: %d\n', N_mics);
fprintf('Microphone spacing (d_y): %.3f m\n', d_y);
fprintf('Source frequency: %d Hz\n', f);
fprintf('Wavelength: %.3f m\n', c_o / f);
fprintf('d_y/wavelength ratio: %.3f\n', d_y / (c_o / f));
fprintf('Source position: (%.3f, %.3f, %.3f) m\n', x_s, y_s, z_s);