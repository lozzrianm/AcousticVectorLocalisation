if ~exist('batch_csv_name', 'var')
    clc; clear all; close all;
else
    clc; close all;
end

% Generate array output using 3D sound field model
% For acoustic vector sensor array
%
% Written by L Marshall 27/11/2025
% Equivalent aperture trial by L Marshall 29/01/2026
% Edited for batch running by L Marshall 20/03/2026

%% CHECK FOR BATCH MODE PARAMETERS %%

batch_mode = exist('batch_csv_name', 'var');
% batch_mode = false;

if batch_mode
    fprintf('\n<strong>BATCH MODE ACTIVE - AVS Signal Generation</strong>\n');
    fprintf('  Frequency: %.0f Hz\n', batch_test_freq);
    fprintf('  Source: (%.3f, %.3f) m\n', batch_test_source_x, batch_test_source_y);
    fprintf('  Output: %s\n', batch_csv_name);
end


%% DEFINE INPUT VARIABLES %%

% Signal Parameters
duration = 10; %sample duration (s)
noise_level = 0.05; %noise amplitude
num_sources = 1; %number of signal sources

% Source positions [x, y, z] (m)
source_positions = [
    -0.4, 0, 0; %source 1
    %4, -2, 0.0; %source 2
];

% Source frequencies (Hz)
source_frequencies = [
    630; %source 1
    %2000; %source 2
];

% Source amplitudes
source_amplitudes = [
    0.01; %source 1
    %0.01; %source 2
];

% Physical Environment Parameters
c_o = 340; %speed of sound (m/s)
rho_o = 1.02; %air density at STP

% Vector Sensor Array Parameters
N_a = 2; %number of independent arrays
N_v = 1; %vector sensors per array
d_y = 0.1; %vector sensor y-axis spacing (m)
delta = 0.04; %MEMS colocation spacing (m)

% Array 1 centre coordinates (m)
x_a1 = 0;
y_a1 = 0;
z_a1 = 0;
A1_coord = [x_a1; y_a1; z_a1];

% Array 2 centre coordinates (m)
x_a2 = 0;
y_a2 = -0.32;
z_a2 = 0;
A2_coord = [x_a2; y_a2; z_a2];

% Sampling Conditions
d_t = 0.0001; %time step (s)
csv_filename = 'generatedsignal_avs.csv'; %default output filename


%% APPLY BATCH MODE PARAMETER OVERRIDES %%

if exist('batch_test_freq', 'var')
    source_frequencies = [batch_test_freq];
end
if exist('batch_test_source_x', 'var') && exist('batch_test_source_y', 'var')
    source_positions = [batch_test_source_x, batch_test_source_y, 0];
    num_sources = 1;
end
if exist('batch_test_delta', 'var')
    delta = batch_test_delta;
end
if exist('batch_csv_name', 'var')
    csv_filename = batch_csv_name;
end


%% INITIALISE ARRAY GEOMETRY %%

% Calculate microphone counts
N_ma = N_v * 4; %MEMS per array (4 per vector sensor)
N_m = N_a * N_ma; %total MEMS in system

% Adjust array centres along y-axis for multi-sensor arrays
if N_v == 1
    y_a1_adjusted = y_a1;
    y_a2_adjusted = y_a2;
else
    y_a1_adjusted = y_a1 - (N_v - 1) * d_y / 2;
    y_a2_adjusted = y_a2 - (N_v - 1) * d_y / 2;
end

% Store array centres as columns [x; y; z]
array_centres = [x_a1, x_a2;
                 y_a1_adjusted, y_a2_adjusted;
                 z_a1, z_a2];


%% DEFINE MICROPHONE POSITIONS %%

mic_positions = zeros(3, N_m);

for array = 1:N_a
    centre = array_centres(:, array);

    for vs = 1:N_v
        idx = (array - 1) * N_ma + (vs - 1) * 4 + (1:4);

        % Centre of this vector sensor along y
        y_c = centre(2) + (vs - 1) * d_y;

        % Four MEMS in square configuration around sensor centre
        mic_positions(:, idx(1)) = [centre(1) - delta/2; y_c - delta/2; centre(3)]; %bottom-left
        mic_positions(:, idx(2)) = [centre(1) + delta/2; y_c - delta/2; centre(3)]; %bottom-right
        mic_positions(:, idx(3)) = [centre(1) + delta/2; y_c + delta/2; centre(3)]; %top-right
        mic_positions(:, idx(4)) = [centre(1) - delta/2; y_c + delta/2; centre(3)]; %top-left
    end
end


%% CREATE TIME VECTOR %%

num_samples = round(duration / d_t);
time = (0:num_samples - 1) * d_t;

fprintf('Generating %d samples over %.1f seconds\n', num_samples, duration);


%% INITIALISE SIGNAL ARRAYS %%

% Dimensions: [time samples, microphones, sources]
p_monopole = zeros(length(time), N_m, num_sources);
p_sinusoidal = zeros(length(time), N_m, num_sources);


%% GENERATE SIGNALS FROM EACH SOURCE %%

fprintf('Generating signals from %d source(s)...\n', num_sources);

for src = 1:num_sources
    fprintf('  Processing source %d/%d (%.1f Hz)...', src, num_sources, source_frequencies(src));

    x_s = source_positions(src, 1);
    y_s = source_positions(src, 2);
    z_s = source_positions(src, 3);
    f_src = source_frequencies(src);
    Q_o = source_amplitudes(src);

    omega = f_src * 2 * pi; %angular frequency

    for k = 1:num_samples
        t1 = time(k);

        for n = 1:N_m
            x_mic = mic_positions(1, n);
            y_mic = mic_positions(2, n);
            z_mic = mic_positions(3, n);

            % Distance from source to microphone
            R_a = sqrt((x_mic - x_s)^2 + (y_mic - y_s)^2 + (z_mic - z_s)^2);

            % Propagation delay
            t = t1 - R_a / c_o;

            Q1_dot = Q_o * omega * cos(omega * t); %monopole derivative
            Q2_dot = Q_o * sin(omega * t); %sinusoidal source

            p_monopole(k, n, src) = (rho_o / (4 * pi * R_a)) * Q1_dot;
            p_sinusoidal(k, n, src) = (rho_o / (4 * pi * R_a)) * Q2_dot;
        end
    end

    fprintf(' Done\n');
end


%% COMBINE SIGNALS FROM ALL SOURCES %%

% Sum contributions across source dimension
p1 = sum(p_monopole, 3);
p2 = sum(p_sinusoidal, 3);

fprintf('Combined signals from all sources\n');


%% ADD SIMULATED NOISE %%

fprintf('Adding noise (level = %.3f)...\n', noise_level);

p1_noisy = p1 + noise_level * randn(size(p1));
p2_noisy = p2 + noise_level * randn(size(p2));


%% NORMALISE AND SAVE TO CSV %%

fprintf('Normalising and saving to CSV...\n');

% Use abs() for robust normalisation across all signal types
generated_monopoles = p1_noisy ./ max(abs(p1_noisy), [], 1);
generated_sinusoidals = p2_noisy ./ max(abs(p2_noisy), [], 1);

results = [time(:), generated_monopoles, generated_sinusoidals];

% Build column headers — monopole channels first, then sinusoidal
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

fprintf('Saved signal data to: %s\n', csv_filename);


%% ARRAY GEOMETRY FIGURE — PUBLICATION STYLE %%
% Single-array and dual-array geometry figures.
% Dashed grey propagation lines run from source to each array centre.
% Microphone labels follow the BL/BR/TR/TL convention used in processing.

if ~batch_mode

    mic_labels = {'M_1', 'M_2', 'M_3', 'M_4'};
    label_offset = 0.01; %label padding (m)
    label_dx = [-1, +1, +1, -1] * label_offset * 2;
    label_dy = [-1, -1, +1, +1] * label_offset * 2;

    % Blue and orange — distinct in greyscale due to luminance difference
    array_colours = [0 0.447 0.741; 0.850 0.325 0.098];


    % SINGLE-ARRAY GEOMETRY FIGURE %
    % 
    % figure('Color', 'w', 'Position', [100 100 560 420]);
    % 
    % ax = axes;
    % set(ax, 'Color', 'w', 'Box', 'on', 'LineWidth', 1.2, ...
    %     'FontName', 'Times New Roman', 'FontSize', 12);
    % hold on; grid on;
    % 
    % % Dashed grey propagation line from source to array 1 centre
    % plot([source_positions(1,1), x_a1], [source_positions(1,2), array_centres(2,1)], ...
    %     '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5, 'HandleVisibility', 'off');
    % 
    % % Array 1 centre as cross marker
    % scatter(x_a1, y_a1, 60, array_colours(1,:), '+', 'LineWidth', 2, ...
    %     'HandleVisibility', 'off');
    % 
    % % Array 1 microphones as filled squares
    % idx_a1 = 1:N_ma;
    % mic_x = mic_positions(1, idx_a1);
    % mic_y = mic_positions(2, idx_a1);
    % 
    % scatter(mic_x, mic_y, 80, array_colours(1,:), 's', 'filled', ...
    %     'DisplayName', 'MEMS Microphone');
    % 
    % % Microphone labels — black text
    % for i = 1:4
    %     text(mic_x(i) + label_dx(i), mic_y(i) + label_dy(i), ...
    %         mic_labels{i}, 'FontName', 'Times New Roman', 'FontSize', 10, ...
    %         'HorizontalAlignment', 'center', 'Color', 'k');
    % end
    % 
    % % Source — bold black X marker
    % scatter(source_positions(1,1), source_positions(1,2), 150, 'kx', ...
    %     'LineWidth', 2.5, 'DisplayName', 'Source');
    % 
    % xlabel('$x$ (m)', 'Interpreter', 'latex', 'FontSize', 14);
    % ylabel('$y$ (m)', 'Interpreter', 'latex', 'FontSize', 14);
    % legend('Location', 'southwest', 'FontName', 'Times New Roman', 'FontSize', 11, 'Box', 'on');
    % set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
    % 
    % % Tight axis limits with minimum padding
    % all_x = [mic_x, source_positions(1,1)];
    % all_y = [mic_y, source_positions(1,2)];
    % x_pad = 0.08 * (max(all_x) - min(all_x));
    % y_pad = 0.08 * (max(all_y) - min(all_y));
    % xlim([min(all_x) - max(x_pad, 0.05), max(all_x) + max(x_pad, 0.05)]);
    % ylim([min(all_y) - max(y_pad, 0.05), max(all_y) + max(y_pad, 0.05)]);
    % axis equal;
    % 
    % exportgraphics(gcf, 'array_geometry_single.pdf', 'ContentType', 'vector');
    % fprintf('Saved single-array geometry figure to: array_geometry_single.pdf\n');


    % DUAL-ARRAY GEOMETRY FIGURE %

    figure('Color', 'w', 'Position', [100 100 560 420]);

    ax = axes;
    set(ax, 'Color', 'w', 'Box', 'on', 'LineWidth', 1.2, ...
        'FontName', 'Times New Roman', 'FontSize', 12);
    hold on; grid on;

    % Dashed grey propagation lines from source to each array centre —
    % one line per array to avoid visual clutter
    plot([source_positions(1,1), x_a1], [source_positions(1,2), array_centres(2,1)], ...
        '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5, 'HandleVisibility', 'off');
    plot([source_positions(1,1), x_a2], [source_positions(1,2), array_centres(2,2)], ...
        '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5, 'HandleVisibility', 'off');

    for array = 1:N_a
        idx_array = (array - 1) * N_ma + (1:N_ma);

        % Array centre as cross marker
        ac = array_centres(1:2, array);
        mic_x = mic_positions(1, idx_array);
        mic_y = mic_positions(2, idx_array);
        scatter(ac(1), ac(2), 60, array_colours(array,:), '+', 'LineWidth', 2, ...
            'HandleVisibility', 'off');

        % Microphones as filled squares
        scatter(mic_x, mic_y, 80, array_colours(array,:), 's', 'filled', ...
            'DisplayName', sprintf('Array %d', array));

        % Microphone labels — black text regardless of array colour
        for i = 1:4
            text(mic_x(i) + label_dx(i), mic_y(i) + label_dy(i), ...
                mic_labels{i}, 'FontName', 'Times New Roman', 'FontSize', 10, ...
                'HorizontalAlignment', 'center', 'Color', 'k');
        end
    end

    % Source — bold black X marker
    for src = 1:num_sources
        scatter(source_positions(src,1), source_positions(src,2), 150, 'kx', ...
            'LineWidth', 2.5, 'DisplayName', 'Source');
    end

    xlabel('$x$ (m)', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('$y$ (m)', 'Interpreter', 'latex', 'FontSize', 14);
    legend('Location', 'southwest', 'FontName', 'Times New Roman', 'FontSize', 11, 'Box', 'on');
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);

    % Tight axis limits with minimum padding
    all_x = [mic_positions(1,:), source_positions(:,1)'];
    all_y = [mic_positions(2,:), source_positions(:,2)'];
    x_pad = 0.08 * (max(all_x) - min(all_x));
    y_pad = 0.08 * (max(all_y) - min(all_y));
    xlim([min(all_x) - max(x_pad, 0.05), max(all_x) + max(x_pad, 0.05)]);
    ylim([min(all_y) - max(y_pad, 0.05), max(all_y) + max(y_pad, 0.05)]);
    axis equal;

    exportgraphics(gcf, 'array_geometry_dual.pdf', 'ContentType', 'vector');
    fprintf('Saved dual-array geometry figure to: array_geometry_dual.pdf\n');

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

array_separation = norm(A2_coord - A1_coord);
fprintf('\nArray 1 centre: (%.3f, %.3f, %.3f) m\n', x_a1, y_a1, z_a1);

if N_a > 1
    fprintf('Array 2 centre: (%.3f, %.3f, %.3f) m\n', x_a2, y_a2, z_a2);
    fprintf('Array separation distance: %.3f m\n', array_separation);
end

if num_sources == 1
    freq_display = source_frequencies(1);
    fprintf('\nSource frequency: %d Hz\n', freq_display);
    fprintf('Wavelength: %.3f m\n', c_o / freq_display);
    fprintf('Delta/wavelength ratio: %.3f\n', delta / (c_o / freq_display));
    fprintf('Array separation/wavelength ratio: %.3f\n', array_separation / (c_o / freq_display));

    if N_v > 1
        fprintf('VS spacing/wavelength ratio: %.3f\n', d_y / (c_o / freq_display));
    end
end

fprintf('\n<strong>Signal Sources: %d</strong>\n', num_sources);
for src = 1:num_sources
    fprintf('Source %d: (%.3f, %.3f, %.3f) m @ %.0f Hz\n', ...
        src, source_positions(src, :), source_frequencies(src));
end