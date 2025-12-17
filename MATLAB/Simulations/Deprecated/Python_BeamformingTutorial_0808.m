clc; clear all; close all;

%L Marshall 05/08/25
%Adapted from the control logic give by Marc Lichtman in PySDR

%% Setting the Inputs %%
Nr = 3; %Number of receiving elements (microphones)
d = 0.5; %Half of wavelength
sample_rate = 1e6; %Sampling rate
N = 10000; %Number of samples
f_tone = 0.02e6; %Tone frequency


time = (0:(N-1))/sample_rate; %Create time vector

tx = exp(2*1i*pi*f_tone*time); %Transmitter signal

n = linspace(0, Nr-1, Nr); %Sensor indices

angles_deg = [0, 20, 45, 60, 90]; %Testing angles


theta_scan = linspace(-pi, pi, 1000); %Scan incoming angles between 180 and -180



%% Sweeping Through Different Angles-of-Attack to Run Sim %%
% Generate inputs signals
combine_sigs = arrayfun(@(theta_deg) ...
    generate_signal(deg2rad(theta_deg), Nr, d, tx), ...
    angles_deg, 'UniformOutput', false);

% Perform beamforming on signals
responses_db = cell2mat(cellfun(@(sig) ...
    beamform_response(sig, Nr, d, theta_scan), ...
    combine_sigs, 'UniformOutput', false));



%% Plot results %%
% Set these ahead of time to simplify functions
num_angles = length(angles_deg);
num_cols = ceil(sqrt(num_angles));
num_rows = ceil(num_angles/num_cols);

% Plot raw input signals
plot_inputsigal(combine_sigs, Nr, angles_deg, num_angles, num_cols, num_rows);

% Plot beamformed arrays
plot_beamforming(theta_scan, responses_db, angles_deg, num_angles, num_cols, num_rows);
plot_polarbeamforming(theta_scan, responses_db, angles_deg, num_angles, num_cols, num_rows);
plot_polarmovie(theta_scan, responses_db, angles_deg, num_angles);



%% Functions %%
% Generate signal for analysis
function combine_sig = generate_signal(theta, Nr, d, tx)
    n = 0:(Nr-1); %Sensor indices
    % Make steering vector then add some noise
    steer_vec = exp(-2*1i*pi*d*n*sin(theta)).'; %Steering vector
    noise = 0.5 * (randn(Nr, length(tx)) + 1i * randn(Nr, length(tx))); %Generate random noise
    combine_sig = steer_vec * tx + noise; %Add noise to signal
end

% Calculate beamforming response
function response_db = beamform_response(combine_sig, Nr, d, theta_scan)
    n = 0:(Nr-1); %Sensor indices
    %Grid values for combos of theta at each mic
    [theta_grid, n_grid] = meshgrid(theta_scan, n); 
    %Weighted matrix
    w_matrix = exp(-2*1i*pi*d .* n_grid .* sin(theta_grid)); 
    beamformed = w_matrix' * combine_sig; %Get get correct matrix dimensions
    power_vals = var(beamformed, 0, 2); %Variance across time
    response_db = 10*log10(power_vals); %Convert to dB (maybe change to 20x?)
    response_db = response_db - max(response_db); %Normalise response
end

% Plot input signals across tested angles
function plot_inputsigal(combine_sigs, Nr, angles_deg, num_angles, num_cols, num_rows)
    figure;
    for k = 1:num_angles
        subplot(num_rows, num_cols, k);
        hold on;
        colors = lines(Nr);
        sig_k = combine_sigs{k};
        for i = 1:Nr
            plot(real(squeeze(sig_k(i,1:200))), 'Color', colors(i,:), ...
                'DisplayName', ['Sensor ' num2str(i)]);
        end
        legend;
        xlabel('Time (samples)');
        ylabel('Amplitude');
        title(['Received Signal at ' num2str(angles_deg(k)) '°']);
        grid on;
    end
end

% Plot beamforming arrays across tested angles
function plot_beamforming(theta_scan, responses_db, angles_deg, ...
    num_angles, num_cols, num_rows)

    figure;
    for k = 1:num_angles
        subplot(num_rows, num_cols, k);
        plot(rad2deg(theta_scan), responses_db(:, k), 'LineWidth', 1);
        title([num2str(angles_deg(k)) '° Angle-of-Attack']);
        xlabel('Theta [Degrees]');
        ylabel('Power [dB]');
        grid on;
        xlim([-180 180]);
    end
    sgtitle('Beamforming Responses for Multiple Angles');
end

%Plot polar beamforming arrays across tested angles
function plot_polarbeamforming(theta_scan, responses_db, angles_deg, ...
    num_angles, num_cols, num_rows)
    figure;

    for k = 1:num_angles
        % Map polar axes to subplot axes then delete the underlying ones
        ax = subplot(num_rows, num_cols, k);
        pax = polaraxes('Parent', ax.Parent);
        pax.Position = ax.Position;
        delete(ax); 

        polarplot(pax, theta_scan, responses_db(:, k), 'LineWidth', 1);
        pax.ThetaZeroLocation = 'top';
        pax.ThetaDir = 'clockwise';
        title(pax, [num2str(angles_deg(k)) '° Angle-of-Attack']);
        grid(pax, 'on');
    end
    sgtitle('Polar Beamforming Responses for Multiple Angles');
end


% Make a movie of the different polar plots
function plot_polarmovie(theta_scan, responses_db, angles_deg, ...
    num_angles)

    fig = figure;
    pax = polaraxes('Parent', fig);
    pax.ThetaZeroLocation = 'top';
    pax.ThetaDir = 'clockwise';
    grid(pax, 'on');

    movie_frames(num_angles) = struct('cdata', [], 'colormap', []);

    for k = 1:num_angles
        cla(pax); %Clear previous plot so that it only captures one subplot per frame
        polarplot(pax, theta_scan, responses_db(:, k), 'LineWidth', 1);
        title(pax, [num2str(angles_deg(k)) '° Angle-of-Attack']);
        drawnow;
        movie_frames(k) = getframe(fig);
    end

    sgtitle('Polar Beamforming Responses Over Time');
    movie(fig, movie_frames, 1, 2);
end

