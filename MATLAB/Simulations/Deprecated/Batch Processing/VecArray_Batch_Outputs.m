clc; clear all; close all;

% Generate array output using 3D sound field model
% For acoustic vector array of two elements each with three subcomponents
% Written by L Marshall 01/08/2025
% Edited by L Marshall 05/10/2025
% Adapted for geometry sweep by ChatGPT 09/10/2025

%% Fixed Source Parameters %%
Q_o = 0.01;
roh_o = 1.02;
c_o = 340; % speed of sound
T = 0.01;
t_0 = 0.02;
d_t = 0.0001; % time step
f = 2000; % frequency (Hz)
omega = f*2*pi;
x_s = 0;
y_s = -0.25;
z_s = 0;

%% Sweep Parameters %%
% Define 10 test configurations for geometry
delta_values = linspace(0.006, 0.012, 4);  % 8-15 mm
d_y_values = linspace(0.04, 0.08, 4);       % 4-8 cm
%d_y_values = linspace(0.04,0.08,1);
config_count = 0;

for dy = d_y_values
    for delt = delta_values
        config_count = config_count + 1;

        fprintf('\n=== Generating configuration %d of 10 ===\n', config_count);
        fprintf('d_y = %.3f m, delta = %.3f m\n', dy, delt);

        %% Acoustic Vector Sensor Array Parameters %%
        x_a = 2; % Array center x-coordinate
        z_a = 0; % Array center z-coordinate

        d_y = dy;
        delta = delt;

        N_v = 3; % Number of vector sensor elements
        N_m = N_v*3; % Total microphones

        %% Define microphone positions %%
        mic_positions = zeros(3, N_m);
        y_centres = [-d_y/2, 0, d_y/2];

        for vs = 1:N_v
            idx = (vs-1)*3+(1:3);
            x_c = x_a;
            y_c = y_centres(vs);
            mic_positions(:,idx(1)) = [x_c; y_c; z_a];
            mic_positions(:,idx(2)) = [x_c + 2*delta; y_c; z_a]; % +x
            mic_positions(:,idx(3)) = [x_c; y_c + delta; z_a]; % +y
        end

        %% Time vector %%
        duration = 60;
        num_samples = round(duration / d_t);
        time = (0:num_samples - 1) * d_t;

        %% Initialise Pressure Arrays %%
        p1 = zeros(length(time), N_m);
        p2 = zeros(length(time), N_m);

        %% Generate signals %%
        for k = 1:num_samples
            t1 = time(k);
            for n = 1:N_m
                x_mic = mic_positions(1, n);
                y_mic = mic_positions(2, n);
                z_mic = mic_positions(3, n);

                R_a = sqrt((x_mic - x_s)^2 + (y_mic - y_s)^2 + (z_mic - z_s)^2);
                t = t1 - R_a/c_o;

                Q1_dot = Q_o*omega*cos(omega*t);
                Q2_dot = Q_o*sin(omega*t);

                p1(k,n) = (roh_o/(4*pi*R_a))*Q1_dot;
                p2(k,n) = (roh_o/(4*pi*R_a))*Q2_dot;
            end
        end

        %% Add Noise %%
        noise_level = 0.05;
        p1_noisy = p1 + noise_level*randn(size(p1));
        p2_noisy = p2 + noise_level*randn(size(p2));

        %% Normalise and Save %%
        generated_monopoles = p1_noisy./max(p1_noisy, [], 1);
        generated_sinusoidals = p2_noisy./max(p2_noisy, [], 1);

        results = [time(:), generated_monopoles, generated_sinusoidals];
        header = {'Time'};

        for vs = 1:N_v
            for mic = 1:3
                header{end+1} = sprintf('AVS%d_Mic%d_Monopole', vs, mic);
            end
        end
        for vs = 1:N_v
            for mic = 1:3
                header{end+1} = sprintf('AVS%d_Mic%d_Sinusoidal', vs, mic);
            end
        end

        results_table = array2table(results, 'VariableNames', header);
        
        % Filename encoding geometry
        filename = sprintf('generatedsignal_avs_dy%.3f_delta%.3f', dy, delt);
        filename = regexprep(filename, '\.', 'p');  % replace decimal dots
        filename = [filename, '.csv'];              % append proper .csv at end
        writetable(results_table, filename);

        %% Display summary for configuration %%
        fprintf('Saved: %s\n', filename);
        fprintf('Delta/wavelength = %.3f\n', delta/(c_o/f));
    end
end

fprintf('\nAll %d configurations generated successfully.\n', config_count)
