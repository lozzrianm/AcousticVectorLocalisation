clc; clear all; close all;

% Compare δ = 0.04 m and δ = 0.06 m sweep results to test whether δ/λ is
% the governing variable for the V₃ Spherical formulation's operating
% envelope.
%
% The hypothesis is that the formulation's behaviour at fixed (r, bearing)
% depends only on the dimensionless ratio δ/λ. If true, two sweeps at
% different δ values plotted against δ/λ should overlay. If they don't
% overlay, there's another length scale at play.
%
% This script:
%   (1) Reads two sweep_results.csv files (one per δ value)
%   (2) Plots each metric vs δ/λ with one line per (bearing, δ) combination
%       — bearing controls colour, δ controls linestyle (solid for first
%       δ, dashed for second) and marker fill (filled vs hollow)
%   (3) Computes RMS overlay error between the two δ curves at matched
%       (bearing) in the δ/λ overlap region, with linear interpolation
%       onto a common δ/λ grid
%   (4) Reports a per-bearing verdict on whether the curves overlay
%
% Inputs: edit RESULTS_FOLDER_DELTA1 and RESULTS_FOLDER_DELTA2 below to
% point to your two sweep result folders. The script picks up the
% sweep_results.csv inside each. It also reads sweep_results.mat for the
% δ value (so you don't have to type it).
%
% Output: figures saved to a new compare_results_<timestamp> folder, plus
% a console summary of the overlay verdict per metric per bearing.
%
% Written by L Marshall 04/05/2026.


%% INPUTS — EDIT THESE %%

RESULTS_FOLDER_DELTA1 = 'SINGLE AVS TEST C 0.04 delta';  % ← edit
RESULTS_FOLDER_DELTA2 = 'SINGLE AVS TEST C 0.06 delta';  % ← edit

% Overlay-quality threshold for the per-bearing verdict. RMS difference
% in log-space (decades) below this is "overlays well", above is "does
% not overlay". Picked at 0.15 ≈ 40% relative variation, which is a
% reasonable bar given the noise of single-realisation MVDR runs.
overlay_threshold_log = 0.15;


%% LOAD DATA %%

fprintf('\n<strong>Loading delta-comparison data</strong>\n');

[d1, label1] = load_sweep(RESULTS_FOLDER_DELTA1);
[d2, label2] = load_sweep(RESULTS_FOLDER_DELTA2);

fprintf('  %s: %d test points\n', label1, height(d1.tab));
fprintf('  %s: %d test points\n', label2, height(d2.tab));

% Sanity checks
if abs(d1.delta - d2.delta) < 1e-6
    error('Both sweeps appear to be at the same delta (%.4f m). Need different deltas to compare.', d1.delta);
end

% Find common bearings
bearings_common = intersect(unique(d1.tab.Bearing_deg), unique(d2.tab.Bearing_deg));
if isempty(bearings_common)
    error('No common bearings between the two sweeps.');
end
fprintf('  Common bearings: %s\n', mat2str(bearings_common'));

% Output folder
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
out_folder = sprintf('compare_results_%s', timestamp);
mkdir(out_folder);
fprintf('  Output folder: %s\n', out_folder);


%% COMPUTE OVERLAP REGION IN δ/λ %%

dl1_range = [min(d1.tab.DeltaLambdaSide), max(d1.tab.DeltaLambdaSide)];
dl2_range = [min(d2.tab.DeltaLambdaSide), max(d2.tab.DeltaLambdaSide)];

dl_overlap_lo = max(dl1_range(1), dl2_range(1));
dl_overlap_hi = min(dl1_range(2), dl2_range(2));

fprintf('\n<strong>δ/λ ranges:</strong>\n');
fprintf('  %s: [%.3f, %.3f]\n', label1, dl1_range);
fprintf('  %s: [%.3f, %.3f]\n', label2, dl2_range);
fprintf('  Overlap: [%.3f, %.3f]\n', dl_overlap_lo, dl_overlap_hi);

if dl_overlap_lo >= dl_overlap_hi
    warning('No δ/λ overlap between sweeps — overlay test cannot be performed. Plots will still be generated.');
    do_overlay_test = false;
else
    do_overlay_test = true;
end


%% PLOT SETUP %%

fn = 'Times New Roman';
fs_ax = 12; fs_lab = 14; fs_leg = 10;
lw = 1.6; ms = 7;

bearing_colours = lines(length(bearings_common));
bearing_markers = {'o', 's', '^', 'd', 'v', '>'};


%% PER-METRIC PLOTS — OVERLAY THE TWO δ VALUES %%

% Each metric: bearing = colour, δ = linestyle.  δ1 = solid+filled
% markers, δ2 = dashed+hollow markers.

metrics_to_plot = struct();
metrics_to_plot(1).field = 'Radial_err_m';
metrics_to_plot(1).ylabel = 'Radial error (m)';
metrics_to_plot(1).log = true;
metrics_to_plot(1).fname = 'overlay_radial_error.png';

metrics_to_plot(2).field = 'Angular_err_deg';
metrics_to_plot(2).ylabel = 'Angular error (deg)';
metrics_to_plot(2).log = true;
metrics_to_plot(2).fname = 'overlay_angular_error.png';

metrics_to_plot(3).field = 'Bracket_3dB_rlambda';
metrics_to_plot(3).ylabel = '3 dB bracket (r/\lambda)';
metrics_to_plot(3).log = true;
metrics_to_plot(3).fname = 'overlay_bracket.png';

metrics_to_plot(4).field = 'FDvel_rel_err';
metrics_to_plot(4).ylabel = '|v_{FD} - v_{true}| / |v_{true}|';
metrics_to_plot(4).log = true;
metrics_to_plot(4).fname = 'overlay_fd_vel_err.png';

metrics_to_plot(5).field = 'FDvel_ang_err_deg';
metrics_to_plot(5).ylabel = 'FD velocity vector angular error (deg)';
metrics_to_plot(5).log = false;
metrics_to_plot(5).fname = 'overlay_fd_ang_err.png';

metrics_to_plot(6).field = 'Beamwidth_3dB_deg';
metrics_to_plot(6).ylabel = '-3 dB beamwidth (deg)';
metrics_to_plot(6).log = false;
metrics_to_plot(6).fname = 'overlay_beamwidth.png';

metrics_to_plot(7).field = 'DI_dB';
metrics_to_plot(7).ylabel = 'Directivity index (dB)';
metrics_to_plot(7).log = false;
metrics_to_plot(7).fname = 'overlay_DI.png';

metrics_to_plot(8).field = 'Num_lobes';
metrics_to_plot(8).ylabel = 'Number of lobes';
metrics_to_plot(8).log = false;
metrics_to_plot(8).fname = 'overlay_num_lobes.png';

% Storage for overlay verdict
verdict_table = struct('metric', {}, 'bearing', {}, 'log_rms', {}, 'verdict', {});


for mi = 1:length(metrics_to_plot)
    m = metrics_to_plot(mi);
    fprintf('\n<strong>Plotting %s</strong>\n', m.field);

    figure('Color', 'w', 'Position', [100 100 700 500]);
    hold on;

    plot_handles = [];
    plot_labels = {};

    for bi = 1:length(bearings_common)
        b_deg = bearings_common(bi);
        col = bearing_colours(bi, :);
        mk = bearing_markers{bi};

        % δ1 series — solid line, filled markers
        mask1 = d1.tab.Bearing_deg == b_deg;
        [dl1, idx1] = sort(d1.tab.DeltaLambdaSide(mask1));
        y1_all = d1.tab.(m.field)(mask1);
        y1 = y1_all(idx1);

        % δ2 series — dashed line, hollow markers
        mask2 = d2.tab.Bearing_deg == b_deg;
        [dl2, idx2] = sort(d2.tab.DeltaLambdaSide(mask2));
        y2_all = d2.tab.(m.field)(mask2);
        y2 = y2_all(idx2);

        if m.log
            y1_plot = max(double(y1), eps);
            y2_plot = max(double(y2), eps);
            h1 = semilogy(dl1, y1_plot, '-', 'Color', col, 'LineWidth', lw, ...
                'Marker', mk, 'MarkerSize', ms, 'MarkerFaceColor', col);
            h2 = semilogy(dl2, y2_plot, '--', 'Color', col, 'LineWidth', lw, ...
                'Marker', mk, 'MarkerSize', ms, 'MarkerFaceColor', 'w');
        else
            h1 = plot(dl1, double(y1), '-', 'Color', col, 'LineWidth', lw, ...
                'Marker', mk, 'MarkerSize', ms, 'MarkerFaceColor', col);
            h2 = plot(dl2, double(y2), '--', 'Color', col, 'LineWidth', lw, ...
                'Marker', mk, 'MarkerSize', ms, 'MarkerFaceColor', 'w');
        end

        plot_handles(end+1) = h1;
        plot_labels{end+1} = sprintf('%.0f°, %s', b_deg, label1);
        plot_handles(end+1) = h2;
        plot_labels{end+1} = sprintf('%.0f°, %s', b_deg, label2);

        % Overlay test in the overlap region
        if do_overlay_test
            in_overlap1 = dl1 >= dl_overlap_lo & dl1 <= dl_overlap_hi;
            in_overlap2 = dl2 >= dl_overlap_lo & dl2 <= dl_overlap_hi;
            if sum(in_overlap1) >= 2 && sum(in_overlap2) >= 2
                % Common δ/λ grid for interpolation
                dl_common = unique([dl1(in_overlap1); dl2(in_overlap2)]);
                dl_common = dl_common(dl_common >= dl_overlap_lo & ...
                                      dl_common <= dl_overlap_hi);
                if length(dl_common) >= 2
                    y1_interp = interp1(dl1, double(y1), dl_common, 'linear');
                    y2_interp = interp1(dl2, double(y2), dl_common, 'linear');

                    % Skip NaN (saturated) points
                    valid = ~isnan(y1_interp) & ~isnan(y2_interp) & ...
                            isfinite(y1_interp) & isfinite(y2_interp);

                    if sum(valid) >= 2
                        if m.log
                            log_rms = sqrt(mean( ...
                                (log10(max(y1_interp(valid), eps)) - ...
                                 log10(max(y2_interp(valid), eps))).^2));
                        else
                            % For linear metrics, normalise by mean
                            mu = mean([y1_interp(valid); y2_interp(valid)]);
                            if abs(mu) > eps
                                log_rms = sqrt(mean( ...
                                    ((y1_interp(valid) - y2_interp(valid)) / mu).^2));
                            else
                                log_rms = NaN;
                            end
                        end

                        if isnan(log_rms)
                            verdict = 'undefined';
                        elseif log_rms < overlay_threshold_log
                            verdict = 'OVERLAYS';
                        else
                            verdict = 'does NOT overlay';
                        end

                        verdict_table(end+1) = struct('metric', m.field, ...
                            'bearing', b_deg, 'log_rms', log_rms, ...
                            'verdict', verdict);

                        fprintf('  Bearing %.0f°: log-RMS = %.3f → %s\n', ...
                            b_deg, log_rms, verdict);
                    end
                end
            end
        end
    end

    % Mark the overlap region
    if do_overlay_test
        yl = ylim;
        patch([dl_overlap_lo dl_overlap_hi dl_overlap_hi dl_overlap_lo], ...
              [yl(1) yl(1) yl(2) yl(2)], [0.9 0.9 0.9], ...
              'EdgeColor', 'none', 'FaceAlpha', 0.2, 'HandleVisibility', 'off');
        % Move patch behind the lines
        chld = get(gca, 'Children');
        set(gca, 'Children', flipud(chld));
        ylim(yl);
    end

    xlabel('\delta / \lambda', 'FontName', fn, 'FontSize', fs_lab);
    ylabel(m.ylabel, 'FontName', fn, 'FontSize', fs_lab);
    legend(plot_handles, plot_labels, 'Location', 'best', ...
        'FontName', fn, 'FontSize', fs_leg, 'Box', 'on');
    set(gca, 'FontName', fn, 'FontSize', fs_ax, 'Box', 'on', 'YGrid', 'on');
    title(sprintf('%s — δ comparison', strrep(m.field, '_', ' ')), ...
        'FontName', fn, 'FontWeight', 'normal');

    exportgraphics(gcf, fullfile(out_folder, m.fname), 'ContentType', 'image');
end


%% OVERLAY VERDICT SUMMARY %%

if do_overlay_test && ~isempty(verdict_table)
    fprintf('\n========================================\n');
    fprintf('<strong>OVERLAY VERDICT SUMMARY</strong>\n');
    fprintf('========================================\n');
    fprintf('Threshold for "overlays well": log-RMS < %.3f\n\n', overlay_threshold_log);

    fprintf('%-30s %-9s %-12s %s\n', 'Metric', 'Bearing', 'log-RMS', 'Verdict');
    fprintf('%s\n', repmat('-', 1, 75));
    for i = 1:length(verdict_table)
        v = verdict_table(i);
        fprintf('%-30s %-9.0f %-12.3f %s\n', v.metric, v.bearing, ...
            v.log_rms, v.verdict);
    end

    % Aggregate verdict per metric (across all bearings)
    fprintf('\n<strong>Aggregate verdict per metric (all bearings):</strong>\n');
    metric_names = unique({verdict_table.metric});
    for mi = 1:length(metric_names)
        match = strcmp({verdict_table.metric}, metric_names{mi});
        rms_vals = [verdict_table(match).log_rms];
        rms_vals = rms_vals(~isnan(rms_vals));
        if isempty(rms_vals)
            fprintf('  %-30s: insufficient overlap data\n', metric_names{mi});
        else
            mean_rms = mean(rms_vals);
            verdict_overall = 'OVERLAYS';
            if mean_rms >= overlay_threshold_log
                verdict_overall = 'does NOT overlay';
            end
            fprintf('  %-30s: mean log-RMS = %.3f → %s\n', ...
                metric_names{mi}, mean_rms, verdict_overall);
        end
    end

    % Save verdict to a CSV
    if ~isempty(verdict_table)
        verdict_t = struct2table(verdict_table);
        writetable(verdict_t, fullfile(out_folder, 'overlay_verdict.csv'));
        fprintf('\nVerdict saved: %s\n', fullfile(out_folder, 'overlay_verdict.csv'));
    end
end


%% INTERPRETATION HINT %%

fprintf('\n<strong>How to interpret these results:</strong>\n');
fprintf('  - If MOST metrics OVERLAY: δ/λ is the governing variable for\n');
fprintf('    those metrics; the formulation behaviour collapses onto a\n');
fprintf('    single curve when plotted against δ/λ.\n');
fprintf('  - If metrics DO NOT OVERLAY: there is another length scale\n');
fprintf('    at play (most likely r/λ, or absolute δ via FD numerics).\n');
fprintf('  - The overlap region is gray-shaded on each plot. Disagreement\n');
fprintf('    inside that region is meaningful; disagreement outside it\n');
fprintf('    just means one sweep didn''t test that δ/λ point.\n');
fprintf('\n  Done. Results in: %s\n\n', out_folder);


%% LOCAL FUNCTIONS %%

function [d, label] = load_sweep(folder)
    csv_path = fullfile(folder, 'sweep_results.csv');
    mat_path = fullfile(folder, 'sweep_results.mat');

    if ~isfile(csv_path)
        error('Cannot find sweep_results.csv in: %s', folder);
    end
    d.tab = readtable(csv_path);

    if isfile(mat_path)
        s = load(mat_path);
        if isfield(s, 'delta_fixed')
            d.delta = s.delta_fixed;
        else
            % Fallback — recover δ from the first row
            % DeltaLambdaSide × Lambda_m = δ
            d.delta = d.tab.DeltaLambdaSide(1) * d.tab.Lambda_m(1);
            warning('No delta_fixed in .mat — recovered %.4f m from first row', d.delta);
        end
    else
        d.delta = d.tab.DeltaLambdaSide(1) * d.tab.Lambda_m(1);
        warning('No .mat file in %s — recovered delta = %.4f m', folder, d.delta);
    end

    label = sprintf('\\delta=%.3fm', d.delta);
end
