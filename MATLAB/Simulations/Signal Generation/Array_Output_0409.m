clc; clear all; close all;

% Array_Output_0409.m

%Generate array output using 3D sound field model
%Adapted from code by J Pan 24/6/2025
%Written by L Marshall 01/08/2025
%Edited by L Marshall 04/09/2025

%% Source parameters %%
Q_o = 0.01;
roh_o = 1.02;
c_o = 340; %speed of sound
T = 0.01;
t_0 = 0.02;
d_t = 0.0001; %time step
f = 1000; %frequency (Hz)
omega = f*2*pi; %angular freq
x_s = 0;
y_s = -0.25;
z_s = 0;

%% Array parameters %%
x_a = 2;
z_a = 0;
d_y = 0.08; %spacing between array elements

N_mics = 11;

%% Time vector for 20 seconds %%
duration = 60; % seconds
num_samples = round(duration / d_t);
time = (0:num_samples - 1) * d_t;

%% Initialise Pressure Arrays %%

p1 = zeros(length(time), N_mics); %Monopole
p2 = zeros(length(time), N_mics); %Sinusoidal

%% time loop %%

for k = 1:num_samples
    %t1 = (k-1)*d_t;
    t1 = time(k);
    %array output
    for n = 1:N_mics %number of microphones = 11
        y_a = (n-6)*d_y;
        R_a = sqrt( (x_a-x_s)^2 +(y_a-y_s)^2 + (z_a-z_s)^2 );
        t = t1 -R_a/c_o;

        % Sinusoidal source derivative signals
        Q1_dot = Q_o*omega*cos(omega*t); %Sinusoidal derivative
        Q2_dot = Q_o*sin(omega*t); %Sinusoidal source

        p1(k,n) = (roh_o/(4*pi*R_a))*Q1_dot;
        p2(k,n) = (roh_o/(4*pi*R_a))*Q2_dot;

    end
end

% plot(time, p1(1,:),'k-','LineWidth', 2);
% plot(time, p2(1,:),'k-','LineWidth', 2);

%% L Marshall additions %% 
%% Add Noise %%

noise_level = 0.05; %adjust amplitude
p1_noisy = p1+noise_level*randn(size(p1));
p2_noisy = p2+noise_level*randn(size(p2));


plot(time, p1_noisy(:,1)/max(p1_noisy(:,1)),'k-',time, p2_noisy(:,1)/max(p2_noisy(:,1)),'r-','LineWidth', 1);
ylabel('Array Output (V)')
xlabel('Time (Seconds)')
legend('Monopole Derivative','Sinusoidal Source')
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14); % optional for report styling
xlim([0 0.01]);
ylim([-1.5 1.5]);


%% Saving Results %%
generated_monopoles = p1_noisy./max(p1_noisy, [], 1);
generated_sinusoidals = p2_noisy./max(p2_noisy, [], 1);
%results = [time(:), generated_monopoles(:), generated_modulateds(:)];
results = [time(:), generated_monopoles, generated_sinusoidals];

% Save data to file
writematrix(results, 'generatedsignal.csv')

% Source and array center positions
x_array = x_a;
y_array = 0;

% Vector from source to array center
dx = x_array - x_s;
dy = y_array - y_s;

% Direction of Arrival (DoA) angle in radians and degrees
theta_rad = atan2(-dy, dx); % Negative dy because y-axis increases upward
theta_deg = rad2deg(theta_rad);

fprintf('Direction of Arrival (DoA): %.2f degrees\n', theta_deg);

%% 2D Plot of Linear Microphone Array Geometry %%
figure;
hold on;
grid on;
axis equal;

% Microphone positions
num_mics = 11;
y_positions = linspace(-0.25, 0.25, num_mics);  % evenly spaced and compact
x_positions = x_a * ones(1, num_mics);          % all mics at x = x_a

% Plot microphones
scatter(x_positions, y_positions, 80, 'b', 'filled');

% Plot source (no label)
scatter(x_s, y_s, 100, 'r', 'filled');

for n = 1:num_mics
    plot([x_s, x_positions(n)], [y_s, y_positions(n)], ...
        '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5);
end

% Axes labels and formatting
xlabel('X (m)');
ylabel('Y (m)');
%title('Layout of Linear Microphone Array');
legend('Microphones', 'Source', 'Location', 'northwest');

% Tidy axis limits for professional spacing
xlim([-0.3, 2.1]);
ylim([-0.4, 0.4]);

set(gca, 'FontName', 'Times New Roman', 'FontSize', 14); % optional for report styling
box on;