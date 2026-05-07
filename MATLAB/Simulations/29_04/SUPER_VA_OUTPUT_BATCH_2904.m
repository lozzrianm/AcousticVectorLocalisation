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
% Corrected signal model and normalisation by L Marshall 29/04/2026
%
% Changes (29/04/2026):
%   - Single monopole-derivative signal per channel (cos+sin pair was
%     redundant at one frequency — both components carry identical
%     spatial information once normalised)
%   - Global power normalisation rather than per-channel peak — preserves
%     the 1/R amplitude weighting that the AVS steering vector explicitly
%     requires (full Green's function, not phase-only)
%   - Source y default offset by -0.005 m for consistency with ULA scripts;
%     AVS cardioid pattern is largely insensitive to this but keeping the
%     same default removes one source of cross-architecture variability
%   - Vectorised the signal generation inner loop


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
noise_level = 0; %noise amplitude (relative to unit-power signal)
num_sources = 1; %number of signal sources

% Source positions [x, y, z] (m) — small y offset for ULA compatibility
source_positions = [
    -0.4, -0.005, 0; %source 1
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
N_a = 2; %number of independent arrays — must match MVDR script setting
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
if exist('batch_test_N_a', 'var')
    N_a = batch_test_N_a;
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
% MEMS layout per vector sensor: BL, BR, TR, TL square at side length delta
% This ordering MUST match the mic_offsets matrix used in the MVDR script

mic_positions = zeros(3, N_m);

for array = 1:N_a
    centre = array_centres(:, array);

    for vs = 1:N_v
        idx = (array - 1) * N_ma + (vs - 1) * 4 + (1:4);

        % Centre of this vector sensor along y
        y_c = centre(2) + (vs - 1) * d_y;

        % Four MEMS in square configuration around sensor centre
        mic_positions(:, idx(1)) = [centre(1) - delta/2; y_c - delta/2; centre(3)]; %BL
        mic_positions(:, idx(2)) = [centre(1) + delta/2; y_c - delta/2; centre(3)]; %BR
        mic_positions(:, idx(3)) = [centre(1) + delta/2; y_c + delta/2; centre(3)]; %TR
        mic_positions(:, idx(4)) = [centre(1) - delta/2; y_c + delta/2; centre(3)]; %TL
    end
end


%% CREATE TIME VECTOR %%

num_samples = round(duration / d_t);
time = (0:num_samples - 1) * d_t;

fprintf('Generating %d samples over %.1f seconds\n', num_samples, duration);


%% GENERATE SIGNAL %%
% Vectorised over time and microphones for each source — single
% monopole-derivative pressure field per channel; the previous cos+sin
% pair was deterministically locked in quadrature and carried no
% additional spatial information at one frequency

fprintf('Generating signal (vectorised) from %d source(s)...\n', num_sources);

p_total = zeros(length(time), N_m);

for src = 1:num_sources
    fprintf('  Processing source %d/%d (%.1f Hz)...', src, num_sources, source_frequencies(src));

    x_s = source_positions(src, 1);
    y_s = source_positions(src, 2);
    z_s = source_positions(src, 3);
    f_src = source_frequencies(src);
    Q_o = source_amplitudes(src);
    omega = f_src * 2 * pi;

    % Source-to-mic distances (1 x N_m)
    R_a = sqrt((mic_positions(1,:) - x_s).^2 + ...
               (mic_positions(2,:) - y_s).^2 + ...
               (mic_positions(3,:) - z_s).^2);

    % Retarded-time argument (num_samples x N_m) — broadcast over mics
    t_ret = time(:) - R_a / c_o;

    % Pressure field — monopole derivative dQ/dt = Q_o * omega * cos(omega*t)
    % Amplitude carries the 1/R Green's function term (paper Eq. 1)
    p_src = (rho_o ./ (4 * pi * R_a)) .* (Q_o * omega * cos(omega * t_ret));

    p_total = p_total + p_src;

    fprintf(' Done\n');
end


%% ADD SIMULATED NOISE %%
% Zero-mean Gaussian noise — uncorrelated across channels, fixed amplitude
% relative to the un-normalised pressure field

fprintf('Adding noise (level = %.3f)...\n', noise_level);

p_noisy = p_total + noise_level * randn(size(p_total));


%% GLOBAL POWER NORMALISATION %%
% Normalise by reference (loudest) channel power — preserves the 1/R
% amplitude weighting between channels rather than flattening it. The AVS
% steering vector includes amplitude terms (full Green's function) so
% per-channel normalisation would invalidate the steering vector match.
% Reference channel is the closest mic to the first source.

x_s_ref = source_positions(1, 1);
y_s_ref = source_positions(1, 2);
z_s_ref = source_positions(1, 3);
R_ref = sqrt((mic_positions(1,:) - x_s_ref).^2 + ...
             (mic_positions(2,:) - y_s_ref).^2 + ...
             (mic_positions(3,:) - z_s_ref).^2);
[~, ref_idx] = min(R_ref);

P_ref = mean(p_noisy(:, ref_idx).^2);
generated_signal = p_noisy / sqrt(P_ref);

fprintf('Reference channel: M%d (R = %.4f m)\n', ref_idx, R_ref(ref_idx));
fprintf('Reference power: %.4e\n', P_ref);


%% SAVE TO CSV %%
% Single-signal-per-channel format: [time, mic_1, mic_2, ..., mic_N_m]

fprintf('Saving to CSV...\n');

results = [time(:), generated_signal];

% Build column headers — one column per microphone
header = {'Time'};
for array = 1:N_a
    for vs = 1:N_v
        for mic = 1:4
            header{end+1} = sprintf('Array%d_AVS%d_Mic%d', array, vs, mic);
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
