%Generate array output using 3D sound field model
%J Pan 24/6/2025
clc
clear all
close all
%Source parameters
Q_o = 0.01;
roh_o = 1.02;
c_o = 340; %speed of sound
T = 0.01;
t_0 = 0.02;
d_t = 0.0001; %time step
omega = 500*2*pi; %angular freq
x_s = 0;
y_s = -0.25;
z_s = 0;

%Array parameters
x_a = 2.0;
z_a = 0;
d_y = 0.08; %spacing between array elements

N_mics = 11;

num_samples = 501;
time = (0:num_samples-1) * d_t; %preallocate how many samples


p1 = zeros(length(time), N_mics);
p2 = zeros(length(time), N_mics);

%time loop

for k = 1:num_samples
    %t1 = (k-1)*d_t;
    t1 = time(k);
    %array output
    for n = 1:N_mics %number of microphones = 11
        y_a = (n-6)*d_y;
        R_a = sqrt( (x_a-x_s)^2 +(y_a-y_s)^2 + (z_a-z_s)^2 );
        t = t1 -R_a/c_o;        
        Q1_dot = - (2*(t - t_0)/T^2)*Q_o*exp(- (t - t_0)^2/T^2);
        Q2_dot = Q_o*exp(-(t - t_0)^2/T^2)*(-(2*(t - t_0)/T^2)*sin(omega*t) + omega*cos(omega*t));
    
        p1(k,n) = (roh_o/(4*pi*R_a))*Q1_dot;
        p2(k,n) = (roh_o/(4*pi*R_a))*Q2_dot;
    end
end

% plot(time, p1(1,:),'k-','LineWidth', 2);
% plot(time, p2(1,:),'k-','LineWidth', 2);

plot(time, p1(:,1)/max(p1(:,1)),'k-',time, p2(:,1)/max(p2(:,1)),'r:','LineWidth', 2);
ylabel('Array Output (V)')
xlabel('Time (Seconds)')
legend('u_1_1','u_2_1')

%% L Marshall additions
generated_monopoles = p1./max(p1, [], 1);
generated_modulateds = p2./max(p2, [], 1);
%results = [time(:), generated_monopoles(:), generated_modulateds(:)];
results = [time(:), generated_monopoles, generated_modulateds];

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
y_positions = ((1:num_mics)-6)*d_y; %Centered around y = 0
x_positions = x_a*ones(1, num_mics); %All mics at x = x_a

% Plot microphones
scatter(x_positions, y_positions, 80, 'b', 'filled');

label_offset_x = 0.05; 
for n = 1:num_mics
    text(x_positions(n) + label_offset_x, y_positions(n), ...
        ['Mic ' num2str(n)], 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
end


% Plot source
scatter(x_s, y_s, 100, 'r', 'filled');
text(x_s - 0.3, y_s + 0.02, 'Source', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');


xlabel('X Position (m)');
ylabel('Y Position (m)');
title('Layout of Linear Microphone Array');
legend('Microphones', 'Source', 'Location', 'northwest');
xlim([-0.35 2.2]);
ylim([-0.6 0.6]);