clc; clear all; close all;

%Generate array output using 3D sound field model
%For acoustic vector array of 1element each with four subcomponents

%Written by L Marshall 01/08/2025
%Edited by L Marshall 05/10/2025
%Edited by L Marshall 22/10/2025

%% Source Parameters %%
Q_o = 0.01;
roh_o = 1.02;
c_o = 340; %speed of sound
T = 0.01;
t_0 = 0.02;
d_t = 0.0001; %time step
f = 2000; %frequency (Hz)
omega = f*2*pi; %angular freq
x_s = 0;
y_s = -0.25;
z_s = 0;

%% Acoustic Vector Sensor Array Parameters %%
% Finite Difference method delta
delta = 0.04; %(m)

% Array element spacing (in y)  
d_y = 0.04; %(m)

% Define geom centre of the array 
x_a = 2.0 + delta/2; % Geometric centre x-coordinate
z_a = 0; % Array centre z-coordinate
y_a = -d_y/2 + delta/2; % Geometric centre y-coordinate

% Four MEMS microphones per sensor
N_v = 1; %Number of vector sensor elements
N_m = N_v*4; %Total MEMS microphones

%% Define microphone positions %%
mic_positions = zeros(3, N_m);

for vs = 1:N_v
    idx = (vs-1)*4+(1:4);
    x_c = x_a; 
    y_c = y_a;
    
    % Place corners around the center
    mic_positions(:,idx(1)) = [x_c - delta/2; y_c - delta/2; z_a]; % p0: bottom-left
    mic_positions(:,idx(2)) = [x_c + delta/2; y_c - delta/2; z_a]; % p1: bottom-right
    mic_positions(:,idx(3)) = [x_c + delta/2; y_c + delta/2; z_a]; % p2: top-right
    mic_positions(:,idx(4)) = [x_c - delta/2; y_c + delta/2; z_a]; % p3: top-left
end

%% Time vector %%
duration = 60; %(seconds)
num_samples = round(duration / d_t);
time = (0:num_samples - 1) * d_t;

%% Initialise Pressure Arrays %%
p1 = zeros(length(time), N_m); %Monopole derivative
p2 = zeros(length(time), N_m); %Sinusoidal source

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

%% Save generated signals to CSV %%
generated_monopoles = p1_noisy./max(p1_noisy, [], 1);
generated_sinusoidals = p2_noisy./max(p2_noisy, [], 1);

results = [time(:), generated_monopoles, generated_sinusoidals];
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

results_table = array2table(results, 'VariableNames', header);
writetable(results_table, 'generatedsignal_avs.csv');

%% Direction of Arrival %%
x_array = x_a;
y_array = 0;
dx = x_array - x_s;
dy = y_array - y_s;
theta_rad = atan2(-dy, dx);
theta_deg = rad2deg(theta_rad);
fprintf('Direction of Arrival (DoA): %.2f degrees\n', theta_deg);

%% Simplified 2D Array Geometry Plot %%
figure;
hold on; grid on; 

marker_style = 'o';
marker_color = [0 0.447 0.741]; % MATLAB default blue
marker_size = 50;

for vs = 1:N_v
    mic_indices = (vs-1)*4 + (1:4);

    scatter(mic_positions(1, mic_indices), ...
            mic_positions(2, mic_indices), ...
            marker_size, marker_color, 'filled');

    text(mic_positions(1, mic_indices(1))-0.05, ...
         mic_positions(2, mic_indices(1))-0.1, ...
         sprintf('MEMS %d', vs), 'FontSize', 10,'FontName', 'Times New Roman');
end

% Source
scatter(x_s, y_s, 50, 'r', 'filled');

xlabel('X (m)');
ylabel('Y (m)');
%title('Acoustic Vector Sensor Array Layout');
legend('Microphones','Source','Location','northwest');
xlim([-0.3, 2.3]);
ylim([-0.4, 0.4]);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
box on;

%% Display Vector Array Configuration %%
fprintf('\n<strong>AVS Array Configuration </strong>\n');
fprintf('Number of AVS elements: %d\n', N_v);
fprintf('Microphones per AVS: 4 (for FD approximation)\n');
fprintf('Total microphones: %d\n', N_m);
fprintf('Finite difference spacing (delta): %.3f m\n', delta);
fprintf('AVS element spacing: %.3f m\n', d_y);
fprintf('Wavelength at %d Hz: %.3f m\n', f, c_o/f);
fprintf('Delta/wavelength ratio: %.3f\n', delta/(c_o/f));