if ~exist('batch_csv_name', 'var')
    clc; clear all; close all;
else
    clc; close all;
end

% Generate array output using 3D sound field model
% Single acoustic vector sensor — 6-DOF steering vector investigation
% Fixed frequency (1 kHz), variable source distance in units of wavelength
%
% Written by L Marshall 27/11/2025
% 6-DOF distance sweep variant by L Marshall 11/04/2026
% Corrected signal model and normalisation by L Marshall 30/04/2026
%
% Changes (30/04/2026) — propagated from SUPER_VA_OUTPUT_BATCH.m:
%   - Single monopole-derivative signal per channel (cos+sin pair was
%     redundant at one frequency — both components carry identical
%     spatial information once normalised, and summing them in the MVDR
%     pre-processing step just rotated the phase)
%   - Global power normalisation rather than per-channel peak —
%     CRITICAL for 6-DOF investigation: V_6 steering vector contains
%     full Green's function (exp(-jkR_n)/R_n with amplitude). Per-channel
%     peak normalisation flattens the 1/R amplitude variation between
%     mics and removes the spatial information that V_6 is meant to
%     exploit. Reference channel is the closest mic to the source
%   - Vectorised the signal generation inner loop


%% CHECK FOR BATCH MODE PARAMETERS %%

batch_mode = exist('batch_csv_name', 'var');

if batch_mode
    fprintf('\n<strong>BATCH MODE ACTIVE - 6DOF Signal Generation</strong>\n');
    fprintf('  Frequency: %.0f Hz\n', batch_test_freq);
    fprintf('  Source distance: %.4f m (%.2f lambda)\n', batch_test_distance, batch_test_distance_lambda);
    fprintf('  Source: (%.3f, %.3f) m\n', batch_test_source_x, batch_test_source_y);
    fprintf('  Output: %s\n', batch_csv_name);
end


%% DEFINE INPUT VARIABLES %%

% Signal Parameters
duration = 10; %sample duration (s)
noise_level = 0.0; %noise amplitude (relative to unit-power signal)
num_sources = 1; %number of signal sources

% Source positions [x, y, z] (m)
source_positions = [
    -0.4, 0, 0; %source 1
];

% Source frequencies (Hz)
source_frequencies = [
    1000; %source 1 — fixed for distance sweep
];

% Source amplitudes
source_amplitudes = [
    0.01; %source 1
];

% Physical Environment Parameters
c_o = 340; %speed of sound (m/s)
rho_o = 1.02; %air density at STP

% Vector Sensor Array Parameters
N_a = 1; %single array only
N_v = 1; %single vector sensor
d_y = 0.1; %vector sensor y-axis spacing (m) — unused for N_v=1
delta = 0.04; %MEMS colocation spacing (m)

% Array centre coordinates (m)
x_a1 = 0;
y_a1 = 0;
z_a1 = 0;
A1_coord = [x_a1; y_a1; z_a1];

% Sampling Conditions
d_t = 0.0001; %time step (s)
csv_filename = 'generatedsignal_6dof.csv'; %default output filename


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

% Store array centre as column [x; y; z]
array_centres = [x_a1; y_a1; z_a1];


%% DEFINE MICROPHONE POSITIONS %%
% MEMS layout per vector sensor: BL, BR, TR, TL square at side length delta
% This ordering MUST match the mic_offsets matrix used in the MVDR script

mic_positions = zeros(3, N_m);

centre = array_centres(:, 1);

for vs = 1:N_v
    idx = (vs - 1) * 4 + (1:4);

    % Centre of this vector sensor along y
    y_c = centre(2);

    % Four MEMS in square configuration around sensor centre
    mic_positions(:, idx(1)) = [centre(1) - delta/2; y_c - delta/2; centre(3)]; %BL
    mic_positions(:, idx(2)) = [centre(1) + delta/2; y_c - delta/2; centre(3)]; %BR
    mic_positions(:, idx(3)) = [centre(1) + delta/2; y_c + delta/2; centre(3)]; %TR
    mic_positions(:, idx(4)) = [centre(1) - delta/2; y_c + delta/2; centre(3)]; %TL
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
    % Amplitude carries the 1/R Green's function term
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
% CRITICAL for 6-DOF investigation: normalise by reference (closest)
% channel power. This preserves the 1/R amplitude weighting between
% channels — V_6 steering vector contains exp(-jkR_n)/R_n with amplitude,
% so per-channel peak normalisation would invalidate V_6 and degrade it
% back to a phase-only formulation.

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
for vs = 1:N_v
    for mic = 1:4
        header{end+1} = sprintf('AVS%d_Mic%d', vs, mic);
    end
end

results_table = array2table(results, 'VariableNames', header);
writetable(results_table, csv_filename);

fprintf('Saved signal data to: %s\n', csv_filename);


%% DISPLAY CONFIGURATION SUMMARY %%

fprintf('\n<strong>Single AVS Configuration — 6-DOF Distance Sweep</strong>\n');
fprintf('Vector sensors: %d\n', N_v);
fprintf('Microphones per vector sensor: 4\n');
fprintf('Total microphones: %d\n', N_m);
fprintf('Finite difference spacing (delta): %.3f m\n', delta);

freq_display = source_frequencies(1);
fprintf('\nSource frequency: %d Hz\n', freq_display);
fprintf('Wavelength: %.3f m\n', c_o / freq_display);
fprintf('Delta/wavelength ratio: %.3f\n', delta / (c_o / freq_display));

fprintf('\n<strong>Signal Sources: %d</strong>\n', num_sources);
for src = 1:num_sources
    fprintf('Source %d: (%.3f, %.3f, %.3f) m @ %.0f Hz\n', ...
        src, source_positions(src, :), source_frequencies(src));
end
