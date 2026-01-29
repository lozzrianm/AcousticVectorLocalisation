clc; clear all; close all;

% L Marshall - Batch Frequency Sweep Testing
% Runs existing signal generation and processing scripts across frequency range

%% BATCH TEST PARAMETERS %%

batch_test_frequencies = linspace(200, 5000, 25);
batch_num_tests = length(batch_test_frequencies);
delta_fixed = 0.04;
c_0 = 340;

% Create results folder
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
results_folder = sprintf('batch_results_%s', timestamp);
mkdir(results_folder);

fprintf('\n<strong>BATCH FREQUENCY SWEEP TEST</strong>\n');
fprintf('Testing %d frequencies from %.0f Hz to %.0f Hz\n', ...
    batch_num_tests, min(batch_test_frequencies), max(batch_test_frequencies));
fprintf('Fixed delta: %.4f m\n', delta_fixed);
fprintf('Results folder: %s\n\n', results_folder);

%% INITIALISE RESULTS STORAGE %%
% Use unique prefix 'batch_res_' to avoid overwrites from child scripts

batch_res_frequency = zeros(batch_num_tests, 1);
batch_res_wavelength = zeros(batch_num_tests, 1);
batch_res_d_over_lambda = zeros(batch_num_tests, 1);
batch_res_radial_error = zeros(batch_num_tests, 1);
batch_res_angular_error = zeros(batch_num_tests, 1);
batch_res_peak_response_db = zeros(batch_num_tests, 1);
batch_res_beamwidth_deg = zeros(batch_num_tests, 1);
batch_res_sidelobe_level_db = zeros(batch_num_tests, 1);
batch_res_directivity_index_db = zeros(batch_num_tests, 1);
batch_res_condition_number = zeros(batch_num_tests, 1);
batch_res_num_lobes_detected = zeros(batch_num_tests, 1);
batch_res_grating_lobe_present = zeros(batch_num_tests, 1); % Use numeric for table compatibility
% Additional metrics from improved beam pattern analysis
batch_res_beamwidth_6dB = zeros(batch_num_tests, 1);
batch_res_beamwidth_10dB = zeros(batch_num_tests, 1);
batch_res_front_to_back_ratio = zeros(batch_num_tests, 1);
batch_res_peak_to_null_ratio = zeros(batch_num_tests, 1);
batch_res_num_significant_lobes = zeros(batch_num_tests, 1);
batch_res_aliasing_severity = zeros(batch_num_tests, 1); % 0=None, 1=Mild, 2=Moderate, 3=Severe


%% RUN BATCH TESTS %%

for batch_test_idx = 1:batch_num_tests
    
    fprintf('\n========================================\n');
    fprintf('<strong>TEST %d/%d: f = %.0f Hz</strong>\n', ...
        batch_test_idx, batch_num_tests, batch_test_frequencies(batch_test_idx));
    fprintf('========================================\n');
    
    batch_current_freq = batch_test_frequencies(batch_test_idx);
    batch_current_lambda = c_0 / batch_current_freq;
    batch_current_d_over_lambda = delta_fixed / batch_current_lambda;
    
    fprintf('Wavelength: %.4f m\n', batch_current_lambda);
    fprintf('d/lambda ratio: %.4f\n', batch_current_d_over_lambda);
    
    % Store parameters
    batch_res_frequency(batch_test_idx) = batch_current_freq;
    batch_res_wavelength(batch_test_idx) = batch_current_lambda;
    batch_res_d_over_lambda(batch_test_idx) = batch_current_d_over_lambda;
    
    %% STEP 1: CLEAR PREVIOUS RUN VARIABLES %%
    
    % Clear variables from previous script runs to prevent stale data
    vars_to_clear = {'all_results', 'beam_metrics', 'grating_detection', ...
                     'mean_condition_number', 'source_frequencies', 'delta', ...
                     'r_vs_arrays', 'response_db_arrays', 'all_cases_results'};
    for v = 1:length(vars_to_clear)
        if exist(vars_to_clear{v}, 'var')
            clear(vars_to_clear{v});
        end
    end
    
    %% STEP 2: RUN SIGNAL GENERATION SCRIPT %%
    
    fprintf('\nGenerating signal at %.0f Hz...\n', batch_current_freq);
    
    % Set parameters for generation script
    batch_test_freq = batch_current_freq;
    batch_test_delta = delta_fixed;
    batch_csv_name = fullfile(results_folder, sprintf('signal_test_%02d_f%04dHz.csv', ...
                              batch_test_idx, round(batch_current_freq)));
    
    % Verify parameters before running
    fprintf('  batch_test_freq = %.0f Hz\n', batch_test_freq);
    fprintf('  batch_test_delta = %.4f m\n', batch_test_delta);
    fprintf('  batch_csv_name = %s\n', batch_csv_name);
    
    % Run generation script
    run('SUPER_VA_Output.m');
    
    fprintf('  Signal generation complete\n');
    
    %% STEP 3: RUN SIGNAL PROCESSING SCRIPT %%
    
    fprintf('\nRunning MVDR beamforming at %.0f Hz...\n', batch_current_freq);
    
    % Set parameters for processing script
    % Note: batch_test_freq and batch_test_delta are already set above
    batch_csv_input = batch_csv_name;
    batch_save_figures = true;
    batch_figure_prefix = fullfile(results_folder, sprintf('test_%02d_f%04dHz', ...
                                   batch_test_idx, round(batch_current_freq)));
    
    % Verify parameters before running
    fprintf('  batch_csv_input = %s\n', batch_csv_input);
    fprintf('  batch_test_freq = %.0f Hz (should match)\n', batch_test_freq);
    
    % Run processing script (use corrected version)
    run('SUPER_VA_MVDR.m');
    
    fprintf('  Processing complete\n');
    
    %% STEP 4: EXTRACT METRICS FROM WORKSPACE %%
    
    % Verify that the frequency was processed correctly
    if exist('source_frequencies', 'var')
        fprintf('\n  Verification: Processing used f = %.0f Hz\n', source_frequencies(1));
        if abs(source_frequencies(1) - batch_current_freq) > 1
            warning('Frequency mismatch detected! Expected %.0f Hz, got %.0f Hz', ...
                    batch_current_freq, source_frequencies(1));
        end
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
        batch_res_beamwidth_6dB(batch_test_idx) = beam_metrics.beamwidth_6dB;
        batch_res_beamwidth_10dB(batch_test_idx) = beam_metrics.beamwidth_10dB;
        batch_res_sidelobe_level_db(batch_test_idx) = beam_metrics.sidelobe_level;
        batch_res_directivity_index_db(batch_test_idx) = beam_metrics.directivity_index;
        batch_res_front_to_back_ratio(batch_test_idx) = beam_metrics.front_to_back_ratio;
        batch_res_peak_to_null_ratio(batch_test_idx) = beam_metrics.peak_to_null_ratio;
    else
        warning('beam_metrics not found for test %d', batch_test_idx);
        batch_res_beamwidth_deg(batch_test_idx) = NaN;
        batch_res_beamwidth_6dB(batch_test_idx) = NaN;
        batch_res_beamwidth_10dB(batch_test_idx) = NaN;
        batch_res_sidelobe_level_db(batch_test_idx) = NaN;
        batch_res_directivity_index_db(batch_test_idx) = NaN;
        batch_res_front_to_back_ratio(batch_test_idx) = NaN;
        batch_res_peak_to_null_ratio(batch_test_idx) = NaN;
    end
    
    % System metrics
    if exist('mean_condition_number', 'var')
        batch_res_condition_number(batch_test_idx) = mean_condition_number;
    else
        warning('mean_condition_number not found for test %d', batch_test_idx);
        batch_res_condition_number(batch_test_idx) = NaN;
    end
    
    % Grating lobe detection 
    if exist('grating_detection', 'var')
        batch_res_num_lobes_detected(batch_test_idx) = grating_detection.num_lobes;
        batch_res_num_significant_lobes(batch_test_idx) = grating_detection.num_significant_lobes;
        batch_res_grating_lobe_present(batch_test_idx) = double(grating_detection.grating_present);
        batch_res_aliasing_severity(batch_test_idx) = grating_detection.severity_score;
    else
        warning('grating_detection not found for test %d', batch_test_idx);
        batch_res_num_lobes_detected(batch_test_idx) = NaN;
        batch_res_num_significant_lobes(batch_test_idx) = NaN;
        batch_res_grating_lobe_present(batch_test_idx) = 0;
        batch_res_aliasing_severity(batch_test_idx) = 0;
    end
    
    % Print test summary
    fprintf('\n<strong>Test %d Summary (f = %.0f Hz):</strong>\n', batch_test_idx, batch_current_freq);
    fprintf('  d/lambda: %.4f\n', batch_res_d_over_lambda(batch_test_idx));
    fprintf('  Radial error: %.4f m\n', batch_res_radial_error(batch_test_idx));
    fprintf('  Beamwidth (-3dB): %.1f deg\n', batch_res_beamwidth_deg(batch_test_idx));
    fprintf('  Beamwidth (-6dB): %.1f deg\n', batch_res_beamwidth_6dB(batch_test_idx));
    fprintf('  Directivity Index: %.2f dB\n', batch_res_directivity_index_db(batch_test_idx));
    fprintf('  Front-to-Back Ratio: %.2f dB\n', batch_res_front_to_back_ratio(batch_test_idx));
    fprintf('  Significant lobes: %d\n', batch_res_num_significant_lobes(batch_test_idx));
    fprintf('  Aliasing severity: %d\n', batch_res_aliasing_severity(batch_test_idx));

    if batch_res_grating_lobe_present(batch_test_idx)
        fprintf('  <strong>WARNING: Spatial aliasing detected!</strong>\n');
    end
    
    % Close all figures to free memory
    close all;
    
end

%% SAVE RESULTS TO CSV %%

fprintf('\n========================================\n');
fprintf('<strong>SAVING RESULTS</strong>\n');
fprintf('========================================\n');

% Verify all arrays have same length before creating table
fprintf('\nVerifying array dimensions:\n');
fprintf('  batch_res_frequency: %d x %d\n', size(batch_res_frequency));
fprintf('  batch_res_wavelength: %d x %d\n', size(batch_res_wavelength));
fprintf('  batch_res_d_over_lambda: %d x %d\n', size(batch_res_d_over_lambda));
fprintf('  batch_res_radial_error: %d x %d\n', size(batch_res_radial_error));
fprintf('  batch_res_angular_error: %d x %d\n', size(batch_res_angular_error));
fprintf('  batch_res_peak_response_db: %d x %d\n', size(batch_res_peak_response_db));
fprintf('  batch_res_beamwidth_deg: %d x %d\n', size(batch_res_beamwidth_deg));
fprintf('  batch_res_sidelobe_level_db: %d x %d\n', size(batch_res_sidelobe_level_db));
fprintf('  batch_res_directivity_index_db: %d x %d\n', size(batch_res_directivity_index_db));
fprintf('  batch_res_condition_number: %d x %d\n', size(batch_res_condition_number));
fprintf('  batch_res_num_lobes_detected: %d x %d\n', size(batch_res_num_lobes_detected));
fprintf('  batch_res_grating_lobe_present: %d x %d\n', size(batch_res_grating_lobe_present));

% Create results table with all metrics as columns, each test as a row
results_table = table(batch_res_frequency, batch_res_wavelength, batch_res_d_over_lambda, ...
                      batch_res_radial_error, batch_res_angular_error, batch_res_peak_response_db, ...
                      batch_res_beamwidth_deg, batch_res_beamwidth_6dB, batch_res_beamwidth_10dB, ...
                      batch_res_sidelobe_level_db, batch_res_directivity_index_db, ...
                      batch_res_front_to_back_ratio, batch_res_peak_to_null_ratio, ...
                      batch_res_condition_number, batch_res_num_lobes_detected, ...
                      batch_res_num_significant_lobes, batch_res_grating_lobe_present, ...
                      batch_res_aliasing_severity, ...
    'VariableNames', {'Frequency_Hz', 'Wavelength_m', 'd_over_lambda', ...
                      'Radial_Error_m', 'Angular_Error_deg', 'Peak_Response_dB', ...
                      'Beamwidth_3dB_deg', 'Beamwidth_6dB_deg', 'Beamwidth_10dB_deg', ...
                      'Sidelobe_Level_dB', 'Directivity_Index_dB', ...
                      'Front_to_Back_Ratio_dB', 'Peak_to_Null_Ratio_dB', ...
                      'Condition_Number', 'Num_Lobes_Total', ...
                      'Num_Significant_Lobes', 'Aliasing_Present', ...
                      'Aliasing_Severity'});

% Display table preview
fprintf('\nResults table preview:\n');
disp(results_table(1:min(5, height(results_table)), :));

% Save to CSV
csv_filename = fullfile(results_folder, 'batch_test_results.csv');
writetable(results_table, csv_filename);
fprintf('\nResults saved to: %s\n', csv_filename);

% Also save as .mat file for easier MATLAB reloading
mat_filename = fullfile(results_folder, 'batch_test_results.mat');
save(mat_filename, 'results_table', 'batch_test_frequencies', 'delta_fixed', 'c_0');
fprintf('Results also saved to: %s\n', mat_filename);

%% GENERATE SUMMARY PLOTS %%

fprintf('\n<strong>GENERATING SUMMARY PLOTS</strong>\n');

fprintf('\n<strong>GENERATING SUMMARY PLOTS</strong>\n');

% Plot 1: Radial error vs d/lambda
figure('Name', 'Radial Error vs d/lambda');
plot(batch_res_d_over_lambda, batch_res_radial_error, 'b-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
xlabel('d/\lambda', 'FontName', 'Times New Roman', 'FontSize', 14);
ylabel('Radial Error (m)', 'FontName', 'Times New Roman', 'FontSize', 14);
title('Localisation Error vs Microphone Spacing Ratio', ...
      'FontName', 'Times New Roman', 'FontSize', 14);
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);

% Add reference lines
hold on;
xline(0.5, '--r', 'd/\lambda = 0.5', 'LineWidth', 1.5, 'LabelOrientation', 'horizontal');
hold off;

saveas(gcf, fullfile(results_folder, 'summary_radial_error.png'));

% Plot 1.5: Angular error vs d/lambda
figure('Name', 'Angular Error vs d/lambda');
plot(batch_res_d_over_lambda, batch_res_angular_error, 'b-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
xlabel('d/\lambda', 'FontName', 'Times New Roman', 'FontSize', 14);
ylabel('Angular Error (deg)', 'FontName', 'Times New Roman', 'FontSize', 14);
title('Angular Error vs Microphone Spacing Ratio', ...
      'FontName', 'Times New Roman', 'FontSize', 14);
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);

% Plot 2: Beamwidth vs frequency
figure('Name', 'Beamwidth vs Frequency');
plot(batch_res_frequency, batch_res_beamwidth_deg, 'g-o', 'LineWidth', 2, 'MarkerFaceColor', 'g');
xlabel('Frequency (Hz)', 'FontName', 'Times New Roman', 'FontSize', 14);
ylabel('Beamwidth (degrees)', 'FontName', 'Times New Roman', 'FontSize', 14);
title('Beamwidth vs Frequency', 'FontName', 'Times New Roman', 'FontSize', 14);
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
saveas(gcf, fullfile(results_folder, 'summary_beamwidth.png'));

% Plot 3: Grating lobe detection vs d/lambda
figure('Name', 'Grating Lobe Analysis');
subplot(2,1,1);
bar(batch_res_d_over_lambda, batch_res_num_significant_lobes);
xlabel('d/\lambda', 'FontName', 'Times New Roman', 'FontSize', 14);
ylabel('Number of Lobes', 'FontName', 'Times New Roman', 'FontSize', 14);
title('Number of Lobes Detected', 'FontName', 'Times New Roman', 'FontSize', 14);
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);

subplot(2,1,2);
stem(batch_res_d_over_lambda, batch_res_grating_lobe_present, 'r', 'LineWidth', 2, 'MarkerFaceColor', 'r');
xlabel('d/\lambda', 'FontName', 'Times New Roman', 'FontSize', 14);
ylabel('Grating Lobe Present', 'FontName', 'Times New Roman', 'FontSize', 14);
title('Grating Lobe Detection', 'FontName', 'Times New Roman', 'FontSize', 14);
ylim([-0.1, 1.1]);
yticks([0, 1]);
yticklabels({'No', 'Yes'});
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);

saveas(gcf, fullfile(results_folder, 'summary_grating_lobes.png'));

% Plot 4: Condition number vs frequency (log scale)
figure('Name', 'Condition Number vs Frequency');
semilogy(batch_res_frequency, batch_res_condition_number, 'm-o', 'LineWidth', 2, 'MarkerFaceColor', 'm');
xlabel('Frequency (Hz)', 'FontName', 'Times New Roman', 'FontSize', 14);
ylabel('Condition Number', 'FontName', 'Times New Roman', 'FontSize', 14);
title('Matrix Condition Number vs Frequency', 'FontName', 'Times New Roman', 'FontSize', 14);
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
saveas(gcf, fullfile(results_folder, 'summary_condition_number.png'));

% Plot 5: Combined performance metrics
figure('Name', 'Combined Performance Metrics', 'Position', [100, 100, 1200, 400]);

subplot(1,3,1);
plot(batch_res_d_over_lambda, batch_res_radial_error, 'b-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
xlabel('d/\lambda', 'FontName', 'Times New Roman', 'FontSize', 12);
ylabel('Radial Error (m)', 'FontName', 'Times New Roman', 'FontSize', 12);
title('Localisation Accuracy', 'FontName', 'Times New Roman', 'FontSize', 12);
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 10);

subplot(1,3,2);
plot(batch_res_d_over_lambda, batch_res_beamwidth_deg, 'g-o', 'LineWidth', 2, 'MarkerFaceColor', 'g');
xlabel('d/\lambda', 'FontName', 'Times New Roman', 'FontSize', 12);
ylabel('Beamwidth (deg)', 'FontName', 'Times New Roman', 'FontSize', 12);
title('Angular Resolution', 'FontName', 'Times New Roman', 'FontSize', 12);
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 10);

subplot(1,3,3);
yyaxis left;
plot(batch_res_d_over_lambda, batch_res_directivity_index_db, 'r-o', 'LineWidth', 2, 'MarkerFaceColor', 'r');
ylabel('Directivity Index (dB)', 'FontName', 'Times New Roman', 'FontSize', 12);
yyaxis right;
plot(batch_res_d_over_lambda, batch_res_sidelobe_level_db, 'c-s', 'LineWidth', 2, 'MarkerFaceColor', 'c');
ylabel('Sidelobe Level (dB)', 'FontName', 'Times New Roman', 'FontSize', 12);
xlabel('d/\lambda', 'FontName', 'Times New Roman', 'FontSize', 12);
title('Beam Quality', 'FontName', 'Times New Roman', 'FontSize', 12);
legend('DI', 'SLL', 'Location', 'best');
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 10);

saveas(gcf, fullfile(results_folder, 'summary_combined.png'));

%% FINAL SUMMARY %%

fprintf('\n========================================\n');
fprintf('<strong>BATCH TEST COMPLETE</strong>\n');
fprintf('========================================\n');

% Find critical d/lambda values
[~, min_error_idx] = min(batch_res_radial_error);
fprintf('\nBest localisation performance:\n');
fprintf('  Frequency: %.0f Hz\n', batch_res_frequency(min_error_idx));
fprintf('  d/lambda: %.4f\n', batch_res_d_over_lambda(min_error_idx));
fprintf('  Radial error: %.4f m\n', batch_res_radial_error(min_error_idx));

% Find grating lobe threshold
grating_idx = find(batch_res_grating_lobe_present, 1, 'first');
if ~isempty(grating_idx)
    fprintf('\nGrating lobe threshold:\n');
    fprintf('  First appears at f = %.0f Hz\n', batch_res_frequency(grating_idx));
    fprintf('  Corresponding d/lambda = %.4f\n', batch_res_d_over_lambda(grating_idx));
else
    fprintf('\nNo grating lobes detected in test range.\n');
end

fprintf('\nAll results saved to: %s\n', results_folder);
fprintf('\n');