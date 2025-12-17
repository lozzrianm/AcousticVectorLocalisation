clc; clear all; close all;

% Batch Processing Script for MVDR Beamforming
% Processes all generated signal CSV files from geometry sweep
% L Marshall - October 2025

%% Configuration %%
mvdr_script_name = 'VecArray_Batch_MVDR.m';
output_dir = 'MVDR_Results';

%% Find all generated CSV files %%
fprintf('=== MVDR Beamforming Batch Processor ===\n\n');

csv_files = dir('generatedsignal_avs_dy*.csv');

if isempty(csv_files)
    error('No CSV files found. Please run the signal generation script first.');
end

if ~isfile(mvdr_script_name)
    error('MVDR script not found: %s', mvdr_script_name);
end

num_files = length(csv_files);
fprintf('Found %d CSV files to process.\n\n', num_files);

%% Create output directory %%
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
    fprintf('Created output directory: %s\n\n', output_dir);
end

%% Process each CSV file %%
% Columns: Filename, Status, Time, MSE, 1D_max_dB, 1D_y, 2D_max_dB, 2D_x, 2D_y, dx_error, dy_error
processing_log = cell(num_files, 11);

for i = 1:num_files
    current_file = csv_files(i).name;
    
    fprintf('========================================\n');
    fprintf('Processing file %d of %d: %s\n', i, num_files, current_file);
    
    % Extract and display geometry parameters
    tokens = regexp(current_file, 'dy(\w+)_delta(\w+)', 'tokens');
    if ~isempty(tokens)
        dy_str = strrep(tokens{1}{1}, 'p', '.');
        delta_str = strrep(tokens{1}{2}, 'p', '.');
        fprintf('Geometry: d_y = %s m, delta = %s m\n', dy_str, delta_str);
    end
    fprintf('========================================\n\n');
    
    try
        % Pass filename to MVDR script via base workspace
        assignin('base', 'input_filename', current_file);

        assignin('base', 'd_y', str2double(dy_str));
        assignin('base', 'delta', str2double(delta_str));
        
        % Run MVDR beamforming
        tic;
        evalin('base', sprintf('run(''%s'')', mvdr_script_name));
        elapsed_time = toc;
        
        fprintf('\n--- Completed in %.2f seconds ---\n', elapsed_time);
        
        % Retrieve results from base workspace
        try
            scan_results = evalin('base', 'scan_results');
            MSE_value = scan_results.MSE;
            scan_1d_max_db = scan_results.scan_1d_max_db;
            scan_1d_y = scan_results.scan_1d_y;
            scan_2d_max_db = scan_results.scan_2d_max_db;
            scan_2d_x = scan_results.scan_2d_x;
            scan_2d_y = scan_results.scan_2d_y;
            dx_error = scan_results.errors_x;
            dy_error = scan_results.errors_y;
        catch
            % Fallback if scan_results not available
            MSE_value = NaN;
            scan_1d_max_db = NaN;
            scan_1d_y = NaN;
            scan_2d_max_db = NaN;
            scan_2d_x = NaN;
            scan_2d_y = NaN;
            dx_error = NaN;
            dy_error = NaN;
        end
        
        % Save all figures
        fig_handles = findall(0, 'Type', 'figure');
        if ~isempty(fig_handles)
            fprintf('Saving %d figures...\n', length(fig_handles));
            base_name = strrep(current_file, '.csv', '');
            for fig_idx = 1:length(fig_handles)
                fig_filename = sprintf('%s/%s_fig%d.png', output_dir, base_name, fig_idx);
                saveas(fig_handles(fig_idx), fig_filename);
            end
        end
        close all;
        
        % Log success with all metrics
        processing_log{i, 1} = current_file;
        processing_log{i, 2} = 'SUCCESS';
        processing_log{i, 3} = elapsed_time;
        processing_log{i, 4} = MSE_value;
        processing_log{i, 5} = scan_1d_max_db;
        processing_log{i, 6} = scan_1d_y;
        processing_log{i, 7} = scan_2d_max_db;
        processing_log{i, 8} = scan_2d_x;
        processing_log{i, 9} = scan_2d_y;
        processing_log{i, 10} = dx_error;
        processing_log{i, 11} = dy_error;
        
        fprintf('Status: SUCCESS\n');
        fprintf('1D Scan: %.2f dB at y = %.3f m\n', scan_1d_max_db, scan_1d_y);
        fprintf('2D Scan: %.2f dB at (x, y) = (%.3f, %.3f) m\n\n', ...
            scan_2d_max_db, scan_2d_x, scan_2d_y);
        
    catch ME
        % Handle errors
        fprintf('\n!!! ERROR processing %s !!!\n', current_file);
        fprintf('Message: %s\n', ME.message);
        if ~isempty(ME.stack)
            fprintf('Location: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
        end
        fprintf('\n');
        
        close all;
        
        % Log failure
        processing_log{i, 1} = current_file;
        processing_log{i, 2} = 'FAILED';
        processing_log{i, 3} = NaN;
        processing_log{i, 4} = NaN;
        processing_log{i, 5} = NaN;
        processing_log{i, 6} = NaN;
        processing_log{i, 7} = NaN;
        processing_log{i, 8} = NaN;
        processing_log{i, 9} = NaN;
        processing_log{i, 10} = NaN;
        processing_log{i, 11} = NaN;
    end
    
    % Clean workspace for next iteration
    clearvars -except i num_files csv_files output_dir processing_log mvdr_script_name
end

%% Generate Summary Report %%
fprintf('\n========================================\n');
fprintf('BATCH PROCESSING COMPLETE\n');
fprintf('========================================\n\n');

% Count successes and failures
successful = sum(strcmp(processing_log(:,2), 'SUCCESS'));
failed = num_files - successful;

fprintf('Summary:\n');
fprintf('  Total files processed: %d\n', num_files);
fprintf('  Successful: %d\n', successful);
fprintf('  Failed: %d\n', failed);

% Calculate statistics for successful runs
if successful > 0
    success_idx = strcmp(processing_log(:,2), 'SUCCESS');
    times = cell2mat(processing_log(success_idx, 3));
    fprintf('\nTiming Statistics:\n');
    fprintf('  Average time: %.2f seconds\n', mean(times));
    fprintf('  Total time: %.2f seconds\n', sum(times));
    fprintf('  Min time: %.2f seconds\n', min(times));
    fprintf('  Max time: %.2f seconds\n', max(times));
    
    % MSE statistics
    mse_values = cell2mat(processing_log(success_idx, 4));
    if any(~isnan(mse_values))
        fprintf('\nMSE Statistics:\n');
        fprintf('  Average MSE: %.4f m²\n', mean(mse_values(~isnan(mse_values))));
        fprintf('  Min MSE: %.4f m²\n', min(mse_values(~isnan(mse_values))));
        fprintf('  Max MSE: %.4f m²\n', max(mse_values(~isnan(mse_values))));
    end
    
    % 1D Scan statistics
    scan_1d_max = cell2mat(processing_log(success_idx, 5));
    scan_1d_y = cell2mat(processing_log(success_idx, 6));
    if any(~isnan(scan_1d_max))
        fprintf('\n1D Scan Statistics:\n');
        fprintf('  Average max response: %.2f dB\n', mean(scan_1d_max(~isnan(scan_1d_max))));
        fprintf('  Average y position: %.3f m\n', mean(scan_1d_y(~isnan(scan_1d_y))));
    end
    
    % 2D Scan statistics
    scan_2d_max = cell2mat(processing_log(success_idx, 7));
    scan_2d_x = cell2mat(processing_log(success_idx, 8));
    scan_2d_y = cell2mat(processing_log(success_idx, 9));
    if any(~isnan(scan_2d_max))
        fprintf('\n2D Scan Statistics:\n');
        fprintf('  Average max response: %.2f dB\n', mean(scan_2d_max(~isnan(scan_2d_max))));
        fprintf('  Average x position: %.3f m\n', mean(scan_2d_x(~isnan(scan_2d_x))));
        fprintf('  Average y position: %.3f m\n', mean(scan_2d_y(~isnan(scan_2d_y))));
    end
    
    % Position error statistics
    dx_errors = cell2mat(processing_log(success_idx, 10));
    dy_errors = cell2mat(processing_log(success_idx, 11));
    if any(~isnan(dx_errors))
        fprintf('\nPosition Error Statistics:\n');
        fprintf('  Average dx error: %.4f m\n', mean(dx_errors(~isnan(dx_errors))));
        fprintf('  Average dy error: %.4f m\n', mean(dy_errors(~isnan(dy_errors))));
        fprintf('  RMS dx error: %.4f m\n', sqrt(mean(dx_errors(~isnan(dx_errors)).^2)));
        fprintf('  RMS dy error: %.4f m\n', sqrt(mean(dy_errors(~isnan(dy_errors)).^2)));
    end
end

% Display detailed results table
fprintf('\n%-40s %-8s %-8s %-8s %-10s %-10s %-10s %-10s %-10s\n', ...
    'Filename', 'Status', 'Time(s)', 'MSE(m²)', '1D_dB', '1D_y(m)', '2D_dB', '2D_x(m)', '2D_y(m)');
fprintf('%s\n', repmat('-', 1, 150));
for i = 1:num_files
    fprintf('%-40s %-8s %-8.2f %-8.4f %-10.2f %-10.3f %-10.2f %-10.3f %-10.3f\n', ...
        processing_log{i,1}, processing_log{i,2}, ...
        processing_log{i,3}, processing_log{i,4}, processing_log{i,5}, ...
        processing_log{i,6}, processing_log{i,7}, processing_log{i,8}, processing_log{i,9});
end

% Save comprehensive processing log
log_table = cell2table(processing_log, ...
    'VariableNames', {'Filename', 'Status', 'ProcessingTime_s', 'MSE_m2', ...
                      'Scan1D_Max_dB', 'Scan1D_y_m', 'Scan2D_Max_dB', ...
                      'Scan2D_x_m', 'Scan2D_y_m', 'Error_dx_m', 'Error_dy_m'});
log_filename = sprintf('%s/processing_log.csv', output_dir);
writetable(log_table, log_filename);

fprintf('\nResults saved to: %s/\n', output_dir);
fprintf('Processing log: %s\n', log_filename);
fprintf('\nBatch processing finished successfully.\n');