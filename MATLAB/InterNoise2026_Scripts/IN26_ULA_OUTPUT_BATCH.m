if ~exist('batch_csv_name', 'var')
    clc; clear all; close all;
else
    clc; close all;
end

% Generate array output using 3D sound field model
% For uniform linear array of N_mics microphones
%
% Adapted from code by J Pan 24/6/2025
% Written by L Marshall 01/08/2025
% Edited by L Marshall 04/09/2025
% Edited for organisation by L Marshall 14/11/2025
% Edited for batch running by L Marshall 29/03/2026

%% CHECK FOR BATCH MODE PARAMETERS %%

batch_mode = exist('batch_csv_name', 'var');
% batch_mode = false;

if batch_mode
    fprintf('\n<strong>BATCH MODE ACTIVE - ULA Signal Generation</strong>\n');
    fprintf('  Frequency: %.0f Hz\n', batch_test_freq);
    fprintf('  Source: (%.3f, %.3f) m\n', batch_test_source_x, batch_test_source_y);
    fprintf('  Output: %s\n', batch_csv_name);
end


%% DEFINE INPUT VARIABLES %%

% Signal Parameters
duration = 10; %sample duration (s)
noise_level = 0.05; %noise amplitude
Q_o = 0.01; %source amplitude
f = 1000; %source frequency (Hz)

% Source Position
x_s = -0.4; %source x coordinate (m)
y_s = 0; %source y coordinate (m)
z_s = 0; %source z coordinate (m)

% Array Parameters — single sub-array (4 mics per recording)
N_mics = 4; %microphones per recording
x_a = 0; %array x coordinate (m)
z_a = 0; %array z coordinate (m)

% Sub-array centre y coordinates
% Sub-array 1 (upper): centred at y_centre_a1
% Sub-array 2 (lower): centred at y_centre_a2
% Combined geometric centre: y_a
d_y = 0.36 / 7; %inter-element spacing (m) — 0.051429 m (was 0.045)
y_centre_a1 = 0.02 - 1.5 * d_y; %sub-array 1 centre (m) — -0.057143 m (was 0)
y_centre_a2 = -0.34 + 1.5 * d_y; %sub-array 2 centre (m) — -0.262857 m (was -0.2025)
y_a = (y_centre_a1 + y_centre_a2) / 2; %combined geometric centre (m) — -0.160000 m (was -0.10125)

% Physical Environment Parameters
c_o = 340; %speed of sound (m/s)
roh_o = 1.02; %air density at STP

% Sampling Conditions
d_t = 0.0001; %time step (s)
csv_filename = 'generatedsignal.csv'; %default output filename


%% APPLY BATCH MODE PARAMETER OVERRIDES %%

if exist('batch_test_freq', 'var')
    f = batch_test_freq;
end
if exist('batch_test_source_x', 'var') && exist('batch_test_source_y', 'var')
    x_s = batch_test_source_x;
    y_s = batch_test_source_y;
end
if exist('batch_csv_name', 'var')
    csv_filename = batch_csv_name;
end


%% INITIALISE SIGNAL CHARACTERISTICS %%

omega = f * 2 * pi; %angular frequency


%% COMPUTE ARRAY GEOMETRY %%

% Mic y positions — descending order (most positive y first) per sub-array,
% matching the physical channel labelling convention used in processing
n_vec = linspace((N_mics-1)/2, -(N_mics-1)/2, N_mics);
y_mics_a1 = y_centre_a1 + d_y * n_vec;
y_mics_a2 = y_centre_a2 + d_y * n_vec;

% Combined 8-element position vectors
y_mics_all = [y_mics_a1, y_mics_a2];
x_mics_all = x_a * ones(1, 2 * N_mics);


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

    for n = 1:N_mics
        % Microphone y position centred on sub-array 1 (y_centre_a1)
        % Signal generation is identical for both sub-arrays — the
        % processing script handles the two recordings separately
        y_mic = y_centre_a1 + (n - (N_mics + 1)/2) * d_y;
        R_a = sqrt((x_a - x_s)^2 + (y_mic - y_s)^2 + (z_a - z_s)^2);
        t = t1 - R_a / c_o;

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

if ~batch_mode

    figure('Name', 'Array Output Time Domain');
    plot(time, p1_noisy(:,1) / max(abs(p1_noisy(:,1))), 'k-', ...
         time, p2_noisy(:,1) / max(abs(p2_noisy(:,1))), 'r-', 'LineWidth', 1);
    ylabel('Array Output (V)', 'FontName', 'Times New Roman', 'FontSize', 14);
    xlabel('Time (Seconds)', 'FontName', 'Times New Roman', 'FontSize', 14);
    legend('Monopole Derivative', 'Sinusoidal Source');
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
    xlim([0 0.01]);
    ylim([-1.5 1.5]);

end


%% NORMALISE AND SAVE TO CSV %%

fprintf('Normalising and saving to CSV...\n');

% Use abs() for robust normalisation across all signal types
generated_monopoles = p1_noisy ./ max(abs(p1_noisy), [], 1);
generated_sinusoidals = p2_noisy ./ max(abs(p2_noisy), [], 1);

results = [time(:), generated_monopoles, generated_sinusoidals];
writematrix(results, csv_filename);

fprintf('Saved signal data to: %s\n', csv_filename);


%% CALCULATE DIRECTION OF ARRIVAL %%

dx = x_a - x_s;
dy = y_a - y_s;
theta_rad = atan2(-dy, dx);
theta_deg_doa = rad2deg(theta_rad);

fprintf('Direction of Arrival (DoA): %.2f degrees\n', theta_deg_doa);


%% ARRAY GEOMETRY FIGURE — PUBLICATION STYLE %%
% Matches the AVS geometry figure style: Times New Roman, latex axis labels,
% filled squares for microphones, bold X for source, dashed grey propagation
% lines, and $M_n$ microphone labels.
%
% The combined 8-mic ULA is constructed from two sequential 4-mic recordings:
%   Sub-array 1 (upper): centred at y_centre_a1
%   Sub-array 2 (lower): centred at y_centre_a2
% Both sub-arrays share x_a = 0 m.

if ~batch_mode

    label_offset_x = 0.02; %mic label offset to the right of marker (m)

    figure('Name', 'ULA Geometry', 'Color', 'w', 'Position', [100 100 560 420]);

    ax = axes;
    set(ax, 'Color', 'w', 'Box', 'on', 'LineWidth', 1.2, ...
        'FontName', 'Times New Roman', 'FontSize', 12);
    hold on; grid on;

    % Dashed grey propagation lines from source to each sub-array centre —
    % two lines only, one per sub-array, to avoid visual clutter
    plot([x_s, x_a], [y_s, y_centre_a1], '--', 'Color', [0.7 0.7 0.7], ...
        'LineWidth', 1.5, 'HandleVisibility', 'off');
    plot([x_s, x_a], [y_s, y_centre_a2], '--', 'Color', [0.7 0.7 0.7], ...
        'LineWidth', 1.5, 'HandleVisibility', 'off');

    % Sub-array 1 microphones — filled squares (legend entry)
    scatter(x_mics_all(1:N_mics), y_mics_a1, 60, [0.850 0.325 0.098], 's', 'filled', ...
        'DisplayName', 'ULA Microphone');

    % Sub-array 2 microphones — same style, legend suppressed
    scatter(x_mics_all(N_mics+1:end), y_mics_a2, 60, [0.850 0.325 0.098], 's', 'filled', ...
        'HandleVisibility', 'off');

    % Source — bold black X marker, matching AVS geometry figure convention
    scatter(x_s, y_s, 120, 'kx', 'LineWidth', 2.5, ...
        'DisplayName', sprintf('Source (%d Hz)', f));

    % Microphone labels — $M_1$ to $M_8$ treating the combined array as a
    % single 8-element ULA, numbered top to bottom (most positive y first)
    for n = 1:2*N_mics
        text(x_mics_all(n) + label_offset_x, y_mics_all(n), ...
             sprintf('$M_%d$', n), ...
             'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 9, ...
             'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
    end

    % Axis labels — latex interpreter to match AVS figure
    xlabel('$x$ (m)', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('$y$ (m)', 'Interpreter', 'latex', 'FontSize', 14);

    % Legend — southwest to avoid overlap with propagation lines
    legend('Location', 'southwest', 'FontName', 'Times New Roman', 'FontSize', 10, 'Box', 'on');

    set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);

    % Tight axis limits — extra right padding for $M_n$ labels,
    % extra left padding for source marker and propagation line origin
    all_x = [x_mics_all, x_s];
    all_y = [y_mics_all, y_s];
    xlim([min(all_x) - 0.08, max(all_x) + 0.12]);
    ylim([min(all_y) - 0.05, max(all_y) + 0.05]);
    axis equal;

    % Export as vector PDF — same format as AVS geometry figure
    exportgraphics(gcf, 'ULA_geom_config.pdf', 'ContentType', 'vector');
    fprintf('Saved geometry figure to: ULA_geom_config.pdf\n');

end


%% DISPLAY CONFIGURATION SUMMARY %%

fprintf('\n<strong>Linear Array Configuration</strong>\n');
fprintf('Number of microphones: %d (2 x %d sub-arrays)\n', 2*N_mics, N_mics);
fprintf('Microphone spacing (d_y): %.5f m\n', d_y);
fprintf('Sub-array 1 centre: y = %.5f m\n', y_centre_a1);
fprintf('Sub-array 2 centre: y = %.5f m\n', y_centre_a2);
fprintf('Combined geometric centre: y = %.5f m\n', y_a);
fprintf('Total aperture: %.5f m\n', max(y_mics_all) - min(y_mics_all));
fprintf('Source frequency: %.0f Hz\n', f);
fprintf('Wavelength: %.4f m\n', c_o / f);
fprintf('d_y/wavelength ratio: %.4f\n', d_y / (c_o / f));
fprintf('Source position: (%.3f, %.3f, %.3f) m\n', x_s, y_s, z_s);