clc; clear all; close all;

% Generate array output using 3D sound field model
% Adapted from code by J Pan 24/6/2025
% Written by L Marshall 01/08/2025
% Edited by L Marshall 04/09/2025
% Edited for organisation by L Marshall 14/11/2025

%% DEFINE INPUT VARIABLES %%

% Signal Parameters
duration = 10; %sample duration (s)
noise_level = 0.05; %noise amplitude
Q_o = 0.01; %source amplitude
f = 2000; %source frequency (Hz)

% Source Position
x_s = -1; %source x coordinate (m)
y_s = 0; %source y coordinate (m)
z_s = 0; %source z coordinate (m)

% Array Parameters
N_mics = 4; %number of microphones
% d_y = 0.02825; %spacing between array elements (m)
d_y = 0.08; %spacing between array elements (m)
x_a = 0; %array centre x coordinate (m)
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
        y_a = (n - (N_mics + 1)/2) * d_y;
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
y_positions = (1:N_mics) - (N_mics + 1)/2;
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


%% PLOT ARRAY GEOMETRY - CONFERENCE FIGURE %%

figure('Name', 'ULA Geometry', 'Color', 'w', ...
    'Position', [100 100 560 420]);

ax = axes;
set(ax, 'Color', 'w', 'Box', 'on', 'LineWidth', 1.2, ...
    'FontName', 'Times New Roman', 'FontSize', 12);
hold on; grid on;

% Microphone positions
y_positions = (N_mics:-1:1) - (N_mics + 1)/2;
y_positions = y_positions * d_y;
x_positions = x_a * ones(1, N_mics);

% Plot dashed lines from source to each microphone
for n = 1:N_mics
    plot([x_s, x_positions(n)], [y_s, y_positions(n)], ...
        '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5, ...
        'HandleVisibility', 'off');
end

% Plot microphones as filled squares
scatter(x_positions, y_positions, 80, 'bs', 'filled', ...
    'DisplayName', 'MEMS Microphone');

% Plot source as filled triangle
scatter(x_s, y_s, 120, 'r^', 'filled', ...
    'MarkerEdgeColor', 'k', ...
    'DisplayName', sprintf('Source (%d Hz)', f));

% Microphone labels
label_offset_x = 0.04;
mic_labels = arrayfun(@(n) sprintf('M_%d', n), 1:N_mics, 'UniformOutput', false);
for n = 1:N_mics
    text(x_positions(n) + label_offset_x, y_positions(n), ...
         ['$' mic_labels{n} '$'], ...
         'Interpreter', 'latex', ...
         'FontName', 'Times New Roman', 'FontSize', 10, ...
         'HorizontalAlignment', 'left', ...
         'VerticalAlignment', 'middle');
end

% Axis labels
xlabel('$x$ (m)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$y$ (m)', 'Interpreter', 'latex', 'FontSize', 14);

% Legend
legend('Location', 'northwest', 'FontName', 'Times New Roman', ...
    'FontSize', 11, 'Box', 'on');

set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);

% Set axis limits BEFORE annotation coordinate conversion
x_pad = 0.15;
y_pad = 0.05;
xlim([x_s - x_pad, x_a + 0.25]);
ylim([min(y_positions) - y_pad, max(y_positions) + y_pad]);
axis equal;

% Force MATLAB to finalise all axis scaling before reading positions
drawnow; pause(0.05);

% Annotation x position (just right of microphone labels)
annotation_x = x_a + 0.13;

% y positions for M1 and M2 only
y1 = y_positions(1); % bottom microphone
y2 = y_positions(2); % one above

% Read axis properties AFTER axis equal and limits are applied
ax_pos = get(ax, 'Position');
xl = xlim(ax);
yl = ylim(ax);

% Convert data coordinates to normalised figure coordinates
x_norm  = ax_pos(1) + ax_pos(3) * (annotation_x - xl(1)) / (xl(2) - xl(1));
y1_norm = ax_pos(2) + ax_pos(4) * (y1 - yl(1)) / (yl(2) - yl(1));
y2_norm = ax_pos(2) + ax_pos(4) * (y2 - yl(1)) / (yl(2) - yl(1));

% Clamp to valid figure range
x_norm  = max(0.01, min(0.99, x_norm));
y1_norm = max(0.01, min(0.99, y1_norm));
y2_norm = max(0.01, min(0.99, y2_norm));
% 
% % Draw double-headed arrow
% annotation('doublearrow', ...
%     [x_norm, x_norm], [y1_norm, y2_norm], ...
%     'LineWidth', 0.9, ...
%     'Head1Style', 'vback2', ...
%     'Head2Style', 'vback2', ...
%     'Head1Width', 5, ...
%     'Head2Width', 5);
% 
% % Delta label
% text(annotation_x + 0.02, (y1 + y2)/2, '$\delta$', ...
%     'Interpreter', 'latex', ...
%     'FontName', 'Times New Roman', 'FontSize', 11, ...
%     'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');