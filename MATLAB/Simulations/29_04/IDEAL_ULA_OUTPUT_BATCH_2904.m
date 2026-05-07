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
% Reverted to ideal 8-mic ULA by L Marshall 17/04/2026
% Corrected signal model and normalisation by L Marshall 29/04/2026
%
% Changes (29/04/2026):
%   - Single monopole-derivative signal per channel (cos+sin pair was
%     redundant at a single frequency — both components carry identical
%     spatial information once normalised)
%   - Global power normalisation rather than per-channel peak — preserves
%     the 1/R amplitude weighting from the paper's Eq. 1
%   - Source y default offset by -0.005 m to break the broadside symmetry
%     that the ULA cannot resolve (matches realistic positional uncertainty
%     of the physical horn-speaker source)


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
noise_level = 0.05; %noise amplitude (relative to unit-power signal)
Q_o = 0.01; %source amplitude
f = 1000; %source frequency (Hz)

% Source Position — small y offset breaks ULA broadside symmetry,
% which is unresolvable by any pressure-only linear array regardless of
% algorithm and reflects realistic positional uncertainty of the source
x_s = -0.4; %source x coordinate (m)
y_s = -0.005; %source y coordinate (m)
z_s = 0; %source z coordinate (m)

% Array Parameters — ideal 8-mic ULA (single recording)
N_mics = 8; %microphones in array
x_a = 0; %array x coordinate (m)
y_a = -0.16; %array centre y coordinate (m)
z_a = 0; %array z coordinate (m)
d_y = 0.36 / 7; %inter-element spacing (m) — 0.051429 m

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

% Mic y positions — descending order (most positive y first),
% matching the physical channel labelling convention used in processing
n_vec = linspace((N_mics-1)/2, -(N_mics-1)/2, N_mics);
y_mics = y_a + d_y * n_vec;
x_mics = x_a * ones(1, N_mics);


%% CREATE TIME VECTOR %%

num_samples = round(duration / d_t);
time = (0:num_samples - 1) * d_t;

fprintf('Generating %d samples over %.1f seconds\n', num_samples, duration);


%% GENERATE SIGNAL %%
% Vectorised over time and microphones — single monopole derivative per
% channel; the previous cos+sin pair was deterministically locked in
% quadrature and carried no additional spatial information at one frequency

fprintf('Generating signal (vectorised)...\n');

% Source-to-mic distances (1 x N_mics)
R_a = sqrt((x_a - x_s).^2 + (y_mics - y_s).^2 + (z_a - z_s).^2);

% Retarded-time argument (num_samples x N_mics) — broadcast over mics
t_ret = time(:) - R_a / c_o;

% Pressure field — monopole derivative dQ/dt = Q_o * omega * cos(omega*t)
% Amplitude carries the 1/R Green's function term from the paper Eq. 1
p = (roh_o ./ (4 * pi * R_a)) .* (Q_o * omega * cos(omega * t_ret));

fprintf('Signal generation complete\n');


%% ADD SIMULATED NOISE %%
% Zero-mean Gaussian noise — uncorrelated across channels, fixed amplitude
% relative to the un-normalised pressure field

fprintf('Adding noise (level = %.3f)...\n', noise_level);

p_noisy = p + noise_level * randn(size(p));


%% GLOBAL POWER NORMALISATION %%
% Normalise by reference (loudest) channel power — preserves the 1/R
% amplitude weighting between channels rather than flattening it.
% Reference channel is the closest mic to the source (highest SNR).

[~, ref_idx] = min(R_a);
P_ref = mean(p_noisy(:, ref_idx).^2);
generated_signal = p_noisy / sqrt(P_ref);

fprintf('Reference channel: M%d (R = %.4f m)\n', ref_idx, R_a(ref_idx));
fprintf('Reference power: %.4e\n', P_ref);


%% PLOT TIME DOMAIN SIGNAL %%

if ~batch_mode

    figure('Name', 'Array Output Time Domain');
    plot(time, generated_signal(:,1), 'k-', 'LineWidth', 1);
    ylabel('Array Output (normalised)', 'FontName', 'Times New Roman', 'FontSize', 14);
    xlabel('Time (Seconds)', 'FontName', 'Times New Roman', 'FontSize', 14);
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
    xlim([0 0.01]);

end


%% SAVE TO CSV %%

fprintf('Saving to CSV...\n');

% Single-signal-per-channel format: [time, mic_1, mic_2, ..., mic_N]
results = [time(:), generated_signal];
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
% line, and $M_n$ microphone labels.

if ~batch_mode

    label_offset_x = 0.02; %mic label offset to the right of marker (m)

    figure('Name', 'ULA Geometry', 'Color', 'w', 'Position', [100 100 560 420]);

    ax = axes;
    set(ax, 'Color', 'w', 'Box', 'on', 'LineWidth', 1.2, ...
        'FontName', 'Times New Roman', 'FontSize', 12);
    hold on; grid on;

    % Dashed grey propagation line from source to array centre
    plot([x_s, x_a], [y_s, y_a], '--', 'Color', [0.7 0.7 0.7], ...
        'LineWidth', 1.5, 'HandleVisibility', 'off');

    % Microphones — filled squares (legend entry)
    scatter(x_mics, y_mics, 60, [0.850 0.325 0.098], 's', 'filled', ...
        'DisplayName', 'ULA Microphone');

    % Source — bold black X marker, matching AVS geometry figure convention
    scatter(x_s, y_s, 120, 'kx', 'LineWidth', 2.5, ...
        'DisplayName', sprintf('Source (%d Hz)', f));

    % Microphone labels — $M_1$ to $M_8$, numbered top to bottom (most positive y first)
    for n = 1:N_mics
        text(x_mics(n) + label_offset_x, y_mics(n), ...
             sprintf('$M_%d$', n), ...
             'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 9, ...
             'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
    end

    % Axis labels — latex interpreter to match AVS figure
    xlabel('$x$ (m)', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('$y$ (m)', 'Interpreter', 'latex', 'FontSize', 14);

    % Legend — southwest to avoid overlap with propagation line
    legend('Location', 'southwest', 'FontName', 'Times New Roman', 'FontSize', 10, 'Box', 'on');

    set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);

    % Tight axis limits — extra right padding for $M_n$ labels,
    % extra left padding for source marker and propagation line origin
    all_x = [x_mics, x_s];
    all_y = [y_mics, y_s];
    xlim([min(all_x) - 0.08, max(all_x) + 0.12]);
    ylim([min(all_y) - 0.05, max(all_y) + 0.05]);
    axis equal;

    % Export as vector PDF — same format as AVS geometry figure
    exportgraphics(gcf, 'ULA_geom_config.pdf', 'ContentType', 'vector');
    fprintf('Saved geometry figure to: ULA_geom_config.pdf\n');

end


%% DISPLAY CONFIGURATION SUMMARY %%

fprintf('\n<strong>Linear Array Configuration</strong>\n');
fprintf('Number of microphones: %d\n', N_mics);
fprintf('Microphone spacing (d_y): %.5f m\n', d_y);
fprintf('Array centre: y = %.5f m\n', y_a);
fprintf('Total aperture: %.5f m\n', max(y_mics) - min(y_mics));
fprintf('Source frequency: %.0f Hz\n', f);
fprintf('Wavelength: %.4f m\n', c_o / f);
fprintf('d_y/wavelength ratio: %.4f\n', d_y / (c_o / f));
fprintf('Source position: (%.3f, %.3f, %.3f) m\n', x_s, y_s, z_s);
