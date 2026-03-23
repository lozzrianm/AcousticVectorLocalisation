clc; clear all; close all;

% L Marshall - Batch Distance Sweep Testing
% Runs existing signal generation and processing scripts at various source distances
% Tests performance at n wavelengths away from the array

%% BATCH TEST PARAMETERS %%

test_frequency = 1500; % Fixed test frequency (Hz)
c_0 = 340;
lambda = c_0 / test_frequency;

% Distance range in wavelengths
min_wavelengths = 0.1;
max_wavelengths = 1;
num_distance_steps = 10;

batch_test_distances_lambda = linspace(min_wavelengths, max_wavelengths, num_distance_steps);
batch_test_distances = batch_test_distances_lambda * lambda;
batch_num_tests = length(batch_test_distances);

delta_fixed = 0.042;
source_angle_deg = 45; % Source angle from +X axis

% Grid resolution - CONSTANT across all tests
grid_resolution = 0.015; % m per grid point

% Processing method selection
% Set to 'standard' for 3-element [p, vx, vy] or '6element' for [p0, p1, p2, p3, vx, vy]
processing_method = 'standard'; % Options: 'standard', '6element'

% Create results folder
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
if strcmp(processing_method, '6element')
    results_folder = sprintf('distance_sweep_6elem_%s', timestamp);
else
    results_folder = sprintf('distance_sweep_%s', timestamp);
end
mkdir(results_folder);

fprintf('\n<strong>BATCH DISTANCE SWEEP TEST</strong>\n');
fprintf('Processing method: %s\n', processing_method);
fprintf('Testing %d distances from %.2f to %.2f wavelengths\n', ...
    batch_num_tests, min_wavelengths, max_wavelengths);
fprintf('Fixed frequency: %.0f Hz (lambda = %.4f m)\n', test_frequency, lambda);
fprintf('Distance range: %.3f to %.3f m\n', min(batch_test_distances), max(batch_test_distances));
fprintf('Fixed delta: %.4f m\n', delta_fixed);
fprintf('Grid resolution: %.4f m\n', grid_resolution);
fprintf('Source angle: %.0f deg\n', source_angle_deg);
fprintf('Results folder: %s\n\n', results_folder);

%% INITIALISE RESULTS STORAGE %%

batch_res_distance_m = zeros(batch_num_tests, 1);
batch_res_distance_lambda = zeros(batch_num_tests, 1);
batch_res_source_x = zeros(batch_num_tests, 1);
batch_res_source_y = zeros(batch_num_tests, 1);
batch_res_radial_error = zeros(batch_num_tests, 1);
batch_res_angular_error = zeros(batch_num_tests, 1);
batch_res_peak_response_db = zeros(batch_num_tests, 1);
batch_res_beamwidth_deg = zeros(batch_num_tests, 1);
batch_res_sidelobe_level_db = zeros(batch_num_tests, 1);
batch_res_directivity_index_db = zeros(batch_num_tests, 1);
batch_res_condition_number = zeros(batch_num_tests, 1);




%% RUN BATCH TESTS %%

for batch_test_idx = 1:batch_num_tests
    
    fprintf('\n========================================\n');
    fprintf('<strong>TEST %d/%d: r = %.3f m (%.2f lambda)</strong>\n', ...
        batch_test_idx, batch_num_tests, batch_test_distances(batch_test_idx), ...
        batch_test_distances_lambda(batch_test_idx));
    fprintf('========================================\n');
    
    batch_current_distance = batch_test_distances(batch_test_idx);
    batch_current_distance_lambda = batch_test_distances_lambda(batch_test_idx);
    
    % Calculate source position
    source_x = batch_current_distance * cosd(source_angle_deg);
    source_y = batch_current_distance * sind(source_angle_deg);
    
    fprintf('Source position: (%.3f, %.3f) m\n', source_x, source_y);
    fprintf('Frequency: %.0f Hz\n', test_frequency);
    
    % Store parameters
    batch_res_distance_m(batch_test_idx) = batch_current_distance;
    batch_res_distance_lambda(batch_test_idx) = batch_current_distance_lambda;
    batch_res_source_x(batch_test_idx) = source_x;
    batch_res_source_y(batch_test_idx) = source_y;
    
    %% CALCULATE SEARCH AREA %%
    
    % Margin around source
    margin = max(1.5 * lambda, 0.5);
    max_extent = batch_current_distance + margin;
    
    % Calculate grid points for constant resolution
    x_span = 2 * max_extent;
    y_span = 2 * max_extent;
    
    current_x_scan_points = round(x_span / grid_resolution) + 1;
    current_y_scan_points = round(y_span / grid_resolution) + 1;
    
    fprintf('Search area: %.3f x %.3f m\n', x_span, y_span);
    fprintf('Grid: %d x %d points (res = %.4f m)\n', ...
        current_x_scan_points, current_y_scan_points, grid_resolution);
    
    %% STEP 1: CLEAR PREVIOUS RUN VARIABLES %%
    
    vars_to_clear = {'all_results', 'beam_metrics', 'grating_detection', ...
                     'mean_condition_number', 'source_frequencies', 'delta', ...
                     'r_vs_arrays', 'response_db_arrays', 'all_cases_results', ...
                     'source_positions'};
    for v = 1:length(vars_to_clear)
        if exist(vars_to_clear{v}, 'var')
            clear(vars_to_clear{v});
        end
    end
    
    %% STEP 2: RUN SIGNAL GENERATION SCRIPT %%
    
    fprintf('\nGenerating signal at r=%.3f m, f=%.0f Hz...\n', batch_current_distance, test_frequency);
    
    % Set parameters for generation script
    batch_test_freq = test_frequency;
    batch_test_delta = delta_fixed;
    batch_test_source_x = source_x;
    batch_test_source_y = source_y;
    batch_csv_name = fullfile(results_folder, sprintf('signal_test_%02d_r%.3fm.csv', ...
                              batch_test_idx, batch_current_distance));
    
    fprintf('  batch_test_freq = %.0f Hz\n', batch_test_freq);
    fprintf('  batch_test_delta = %.4f m\n', batch_test_delta);
    fprintf('  batch_test_source_x = %.3f m\n', batch_test_source_x);
    fprintf('  batch_test_source_y = %.3f m\n', batch_test_source_y);
    fprintf('  batch_csv_name = %s\n', batch_csv_name);
    
    % Run generation script
    run('SUPER_VA_Output_Distance.m');
    
    fprintf('  Signal generation complete\n');
    
    %% STEP 3: RUN SIGNAL PROCESSING SCRIPT %%
    
    fprintf('\nRunning MVDR beamforming (%s)...\n', processing_method);
    
    % Set parameters for processing script
    batch_csv_input = batch_csv_name;
    batch_save_figures = true;
    batch_figure_prefix = fullfile(results_folder, sprintf('test_%02d_r%.3fm', ...
                                   batch_test_idx, batch_current_distance));
    
    % Override grid parameters
    batch_x_scan_points = current_x_scan_points;
    batch_y_scan_points = current_y_scan_points;
    batch_x_margin = max_extent;
    batch_y_margin = max_extent;
    
    fprintf('  batch_csv_input = %s\n', batch_csv_input);
    fprintf('  batch_test_freq = %.0f Hz\n', batch_test_freq);
    fprintf('  Grid: %d x %d points, margin = %.3f m\n', ...
        batch_x_scan_points, batch_y_scan_points, batch_x_margin);
    
    % Run appropriate processing script
    if strcmp(processing_method, '6element')
        run('SUPER_VA_MVDR_Distance_6element.m');
    else
        run('SUPER_VA_MVDR_Distance.m');
    end
    
    fprintf('  Processing complete\n');
    
    %% STEP 4: EXTRACT METRICS FROM WORKSPACE %%
    
    % Verify processing
    if exist('source_frequencies', 'var')
        fprintf('\n  Verification: f = %.0f Hz, source at (%.3f, %.3f)\n', ...
                source_frequencies(1), source_positions(1,1), source_positions(1,2));
    end
    
    % Localisation metrics
    if exist('all_results', 'var')
        batch_res_radial_error(batch_test_idx) = all_results.radial_error(1);
        batch_res_angular_error(batch_test_idx) = all_results.angular_error(1);
        batch_res_peak_response_db(batch_test_idx) = all_results.peak_response(1);
    else
        warning('all_results not found for test %d', batch_test_idx);
        batch_res_radial_error(batch_test_idx) = NaN;
        batch_res_angular_error(batch_test_idx) = NaN;
        batch_res_peak_response_db(batch_test_idx) = NaN;
    end
    
    % Beam pattern metrics
    if exist('beam_metrics', 'var')
        batch_res_beamwidth_deg(batch_test_idx) = beam_metrics.beamwidth_3dB;
        batch_res_sidelobe_level_db(batch_test_idx) = beam_metrics.sidelobe_level;
        batch_res_directivity_index_db(batch_test_idx) = beam_metrics.directivity_index;
    else
        warning('beam_metrics not found for test %d', batch_test_idx);
        batch_res_beamwidth_deg(batch_test_idx) = NaN;
        batch_res_sidelobe_level_db(batch_test_idx) = NaN;
        batch_res_directivity_index_db(batch_test_idx) = NaN;
    end
    
    % System metrics
    if exist('mean_condition_number', 'var')
        batch_res_condition_number(batch_test_idx) = mean_condition_number;
    else
        warning('mean_condition_number not found for test %d', batch_test_idx);
        batch_res_condition_number(batch_test_idx) = NaN;
    end
    
    % Print test summary
    fprintf('\n<strong>Test %d Summary (r = %.3f m = %.2f lambda):</strong>\n', ...
        batch_test_idx, batch_current_distance, batch_current_distance_lambda);
    fprintf('  Source: (%.3f, %.3f) m\n', source_x, source_y);
    fprintf('  Radial error: %.4f m (%.2f%% of distance)\n', ...
        batch_res_radial_error(batch_test_idx), ...
        100 * batch_res_radial_error(batch_test_idx) / batch_current_distance);
    fprintf('  Angular error: %.3f deg\n', batch_res_angular_error(batch_test_idx));
    fprintf('  Beamwidth (-3dB): %.1f deg\n', batch_res_beamwidth_deg(batch_test_idx));
    fprintf('  Directivity Index: %.2f dB\n', batch_res_directivity_index_db(batch_test_idx));
    
    % Close all figures to free memory
    close all;
    
end

%% SAVE RESULTS TO CSV %%

fprintf('\n========================================\n');
fprintf('<strong>SAVING RESULTS</strong>\n');
fprintf('========================================\n');

% Create results table
results_table = table(batch_res_distance_m, batch_res_distance_lambda, ...
                      batch_res_source_x, batch_res_source_y, ...
                      batch_res_radial_error, batch_res_angular_error, batch_res_peak_response_db, ...
                      batch_res_beamwidth_deg, batch_res_sidelobe_level_db, batch_res_directivity_index_db, ...
                      batch_res_condition_number, ...
    'VariableNames', {'Distance_m', 'Distance_Wavelengths', ...
                      'Source_X_m', 'Source_Y_m', ...
                      'Radial_Error_m', 'Angular_Error_deg', 'Peak_Response_dB', ...
                      'Beamwidth_3dB_deg', 'Sidelobe_Level_dB', 'Directivity_Index_dB', ...
                      'Condition_Number'});

% Display table preview
fprintf('\nResults table preview:\n');
disp(results_table(1:min(5, height(results_table)), :));

% Save to CSV
csv_filename = fullfile(results_folder, 'distance_sweep_results.csv');
writetable(results_table, csv_filename);
fprintf('\nResults saved to: %s\n', csv_filename);

% Also save as .mat file
mat_filename = fullfile(results_folder, 'distance_sweep_results.mat');
save(mat_filename, 'results_table', 'batch_test_distances', 'batch_test_distances_lambda', ...
     'test_frequency', 'lambda', 'delta_fixed', 'c_0', 'source_angle_deg', 'grid_resolution');
fprintf('Results also saved to: %s\n', mat_filename);

%% GENERATE SUMMARY PLOTS %%

fprintf('\n<strong>GENERATING SUMMARY PLOTS</strong>\n');

% Plot 1: Radial error vs distance
figure('Name', 'Radial Error vs Distance');
subplot(2,1,1);
plot(batch_res_distance_m, batch_res_radial_error, 'b-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
xlabel('Distance (m)', 'FontName', 'Times New Roman', 'FontSize', 14);
ylabel('Radial Error (m)', 'FontName', 'Times New Roman', 'FontSize', 14);
title('Radial Error vs Source Distance', 'FontName', 'Times New Roman', 'FontSize', 14);
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);

subplot(2,1,2);
plot(batch_res_distance_lambda, batch_res_radial_error, 'r-o', 'LineWidth', 2, 'MarkerFaceColor', 'r');
xlabel('Distance (wavelengths)', 'FontName', 'Times New Roman', 'FontSize', 14);
ylabel('Radial Error (m)', 'FontName', 'Times New Roman', 'FontSize', 14);
title('Radial Error vs Source Distance (Wavelengths)', 'FontName', 'Times New Roman', 'FontSize', 14);
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);

saveas(gcf, fullfile(results_folder, 'summary_radial_error.png'));

% Plot 2: Relative error
figure('Name', 'Relative Radial Error');
relative_error_pct = 100 * batch_res_radial_error ./ batch_res_distance_m;
plot(batch_res_distance_lambda, relative_error_pct, 'g-o', 'LineWidth', 2, 'MarkerFaceColor', 'g');
xlabel('Distance (wavelengths)', 'FontName', 'Times New Roman', 'FontSize', 14);
ylabel('Relative Error (%)', 'FontName', 'Times New Roman', 'FontSize', 14);
title('Relative Radial Error vs Source Distance', 'FontName', 'Times New Roman', 'FontSize', 14);
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
saveas(gcf, fullfile(results_folder, 'summary_relative_error.png'));

% Plot 3: Angular error
figure('Name', 'Angular Error vs Distance');
plot(batch_res_distance_lambda, batch_res_angular_error, 'm-o', 'LineWidth', 2, 'MarkerFaceColor', 'm');
xlabel('Distance (wavelengths)', 'FontName', 'Times New Roman', 'FontSize', 14);
ylabel('Angular Error (deg)', 'FontName', 'Times New Roman', 'FontSize', 14);
title('Angular Error vs Source Distance', 'FontName', 'Times New Roman', 'FontSize', 14);
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
saveas(gcf, fullfile(results_folder, 'summary_angular_error.png'));

% Plot 4: Combined metrics
figure('Name', 'Combined Performance Metrics', 'Position', [100, 100, 1200, 800]);

subplot(2,2,1);
plot(batch_res_distance_lambda, batch_res_radial_error, 'b-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
xlabel('Distance (\lambda)', 'FontName', 'Times New Roman', 'FontSize', 12);
ylabel('Radial Error (m)', 'FontName', 'Times New Roman', 'FontSize', 12);
title('Localisation Accuracy', 'FontName', 'Times New Roman', 'FontSize', 12);
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 10);

subplot(2,2,2);
plot(batch_res_distance_lambda, batch_res_angular_error, 'm-o', 'LineWidth', 2, 'MarkerFaceColor', 'm');
xlabel('Distance (\lambda)', 'FontName', 'Times New Roman', 'FontSize', 12);
ylabel('Angular Error (deg)', 'FontName', 'Times New Roman', 'FontSize', 12);
title('Angular Accuracy', 'FontName', 'Times New Roman', 'FontSize', 12);
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 10);

subplot(2,2,3);
plot(batch_res_distance_lambda, batch_res_beamwidth_deg, 'g-o', 'LineWidth', 2, 'MarkerFaceColor', 'g');
xlabel('Distance (\lambda)', 'FontName', 'Times New Roman', 'FontSize', 12);
ylabel('Beamwidth (deg)', 'FontName', 'Times New Roman', 'FontSize', 12);
title('Angular Resolution', 'FontName', 'Times New Roman', 'FontSize', 12);
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 10);

subplot(2,2,4);
plot(batch_res_distance_lambda, batch_res_directivity_index_db, 'r-o', 'LineWidth', 2, 'MarkerFaceColor', 'r');
xlabel('Distance (\lambda)', 'FontName', 'Times New Roman', 'FontSize', 12);
ylabel('Directivity Index (dB)', 'FontName', 'Times New Roman', 'FontSize', 12);
title('Directivity', 'FontName', 'Times New Roman', 'FontSize', 12);
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 10);

saveas(gcf, fullfile(results_folder, 'summary_combined.png'));

%% FINAL SUMMARY %%

fprintf('\n========================================\n');
fprintf('<strong>DISTANCE SWEEP TEST COMPLETE</strong>\n');
fprintf('========================================\n');

% Find best performance
[min_error, min_idx] = min(batch_res_radial_error);
fprintf('\nBest localisation performance:\n');
fprintf('  Distance: %.3f m (%.2f lambda)\n', ...
    batch_res_distance_m(min_idx), batch_res_distance_lambda(min_idx));
fprintf('  Radial error: %.4f m\n', min_error);

fprintf('\nRadial Error Statistics:\n');
fprintf('  Mean: %.4f m\n', mean(batch_res_radial_error));
fprintf('  Std Dev: %.4f m\n', std(batch_res_radial_error));

fprintf('\nAll results saved to: %s\n', results_folder);
fprintf('\n');