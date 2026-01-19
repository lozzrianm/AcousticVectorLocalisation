clc; clear all; close all;

% L Marshall - Batch Frequency Sweep Testing
% Runs existing signal generation and processing scripts across frequency range

%% BATCH TEST PARAMETERS %%

test_frequencies = linspace(600, 4600, 21);
num_tests = length(test_frequencies);
delta_fixed = 0.04;
c_0 = 340;

% Create results folder
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
results_folder = sprintf('batch_results_%s', timestamp);
mkdir(results_folder);

fprintf('\n<strong>BATCH FREQUENCY SWEEP TEST</strong>\n');
fprintf('Testing %d frequencies from %.0f Hz to %.0f Hz\n', ...
    num_tests, min(test_frequencies), max(test_frequencies));
fprintf('Results folder: %s\n\n', results_folder);

%% INITIALISE RESULTS STORAGE %%

results = struct();
results.frequency = zeros(num_tests, 1);
results.wavelength = zeros(num_tests, 1);
results.d_over_lambda = zeros(num_tests, 1);
results.radial_error = zeros(num_tests, 1);
results.angular_error = zeros(num_tests, 1);
results.peak_response_db = zeros(num_tests, 1);
results.beamwidth_deg = zeros(num_tests, 1);
results.sidelobe_level_db = zeros(num_tests, 1);
results.directivity_index_db = zeros(num_tests, 1);
results.condition_number = zeros(num_tests, 1);
results.num_lobes_detected = zeros(num_tests, 1);
results.grating_lobe_present = false(num_tests, 1);

%% RUN BATCH TESTS %%

for test_idx = 1:num_tests
    
    fprintf('\n========================================\n');
    fprintf('<strong>TEST %d/%d: f = %.0f Hz</strong>\n', ...
        test_idx, num_tests, test_frequencies(test_idx));
    fprintf('========================================\n');
    
    current_freq = test_frequencies(test_idx);
    lambda = c_0 / current_freq;
    d_over_lambda = delta_fixed / lambda;
    
    fprintf('d/lambda ratio: %.4f\n', d_over_lambda);
    
    % Store parameters
    results.frequency(test_idx) = current_freq;
    results.wavelength(test_idx) = lambda;
    results.d_over_lambda(test_idx) = d_over_lambda;
    
    %% STEP 1: RUN SIGNAL GENERATION SCRIPT %%
    
    fprintf('\nGenerating signal...\n');
    
    % Set parameters for generation script
    batch_test_freq = current_freq;
    batch_test_delta = delta_fixed;
    batch_csv_name = fullfile(results_folder, sprintf('signal_test_%02d.csv', test_idx));
    
    % Run generation script
    run('SUPER_VA_Output.m');
    
    fprintf('  Signal generation complete\n');
    
    %% STEP 2: RUN SIGNAL PROCESSING SCRIPT %%
    
    fprintf('\nRunning MVDR beamforming...\n');
    
    % Set parameters for processing script
    batch_csv_input = batch_csv_name;
    batch_save_figures = true;
    batch_figure_prefix = fullfile(results_folder, sprintf('test_%02d', test_idx));
    
    % Run processing script
    run('SUPER_VA_MVDR.m');
    
    fprintf('  Processing complete\n');
    
    %% STEP 3: EXTRACT METRICS FROM WORKSPACE %%
    
    % Localisation metrics (from all_results structure)
    results.radial_error(test_idx) = all_results.radial_error(1);
    results.angular_error(test_idx) = all_results.angular_error(1);
    results.peak_response(test_idx) = all_results.peak_response(1);
    
    % Beam pattern metrics (from beam_metrics structure)
    results.beamwidth_deg(test_idx) = beam_metrics.beamwidth;
    results.sidelobe_level_db(test_idx) = beam_metrics.sidelobe_level;
    results.directivity_index_db(test_idx) = beam_metrics.directivity_index;
    
    % System metrics
    results.condition_number(test_idx) = mean_condition_number;
    results.num_lobes_detected(test_idx) = grating_detection.num_lobes;
    results.grating_lobe_present(test_idx) = grating_detection.grating_present;
    
    fprintf('\n<strong>Test Summary:</strong>\n');
    fprintf('  Radial error: %.4f m\n', results.radial_error(test_idx));
    fprintf('  Condition number: %.2e\n', results.condition_number(test_idx));
    fprintf('  Num lobes: %d\n', results.num_lobes_detected(test_idx));
    
    if results.grating_lobe_present(test_idx)
        fprintf('  <strong>WARNING: Grating lobes detected!</strong>\n');
    end
    
    % Close all figures to free memory
    close all;
    
end

%% SAVE RESULTS %%

fprintf('\n========================================\n');
fprintf('<strong>SAVING RESULTS</strong>\n');
fprintf('========================================\n');

results_table = struct2table(results);
csv_filename = fullfile(results_folder, 'batch_test_results.csv');
writetable(results_table, csv_filename);

fprintf('Results saved to: %s\n', csv_filename);
fprintf('\n<strong>BATCH TEST COMPLETE</strong>\n\n');