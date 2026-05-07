clc; clear all; close all;

% Replot delta x distance sweep figures from results CSV
%
% Standalone script: reads delta_dist_results.csv produced by
% BATCH_DELTA_DIST_SWEEP_2D.m and regenerates the publication
% figures with thesis-quality formatting. Does NOT re-run any
% MVDR processing — purely a postprocessing / styling pass on
% existing results.
%
% Use this script to iterate on figure style (fonts, axis labels,
% colour maps, sizes) without re-running the full ~50-minute
% 2D sweep. Designed for thesis main-body and appendix output.
%
% Outputs (in same folder as the CSV):
%   crossover_map_thesis.png/pdf/fig
%   saturation_boundary_thesis.png/pdf/fig
%   velocity_error_map_thesis.png/pdf/fig
%   bracket_grid_thesis.png/pdf/fig
%
% Written by L Marshall 03/05/2026


%% USER INPUTS %%

% Path to the results CSV from BATCH_DELTA_DIST_SWEEP_2D.m
csv_path = 'delta_dist_results.csv';

% Test parameters (must match the run that produced the CSV)
test_freq = 1000; %Hz
c_0 = 340; %m/s
lambda = c_0 / test_freq;

% Operational boundaries (from §3.4 of InterNoise 2026 paper)
velocity_lower_bound = 0.07;
velocity_upper_bound = 0.33;
alias_threshold = 0.5;


%% LOAD AND ARRANGE DATA %%

fprintf('Loading results from: %s\n', csv_path);

T = readtable(csv_path);

delta_lambdas = unique(T.Delta_Lambda);
distance_lambdas = unique(T.Distance_Lambda);
num_delta = length(delta_lambdas);
num_dist = length(distance_lambdas);
num_methods = 4;

fprintf('  delta/lambda values (%d): ', num_delta);
fprintf('%.2f  ', delta_lambdas); fprintf('\n');
fprintf('  r/lambda values (%d): ', num_dist);
fprintf('%.2f  ', distance_lambdas); fprintf('\n');

% Build [num_delta x num_dist] tensors keyed by (Delta_Lambda, Distance_Lambda)
bracket_planar    = build_tensor(T, 'Bracket_Planar_m',    delta_lambdas, distance_lambdas);
bracket_spherical = build_tensor(T, 'Bracket_Spherical_m', delta_lambdas, distance_lambdas);
bracket_modified  = build_tensor(T, 'Bracket_Modified_m',  delta_lambdas, distance_lambdas);
bracket_presonly  = build_tensor(T, 'Bracket_PresOnly_m',  delta_lambdas, distance_lambdas);
sat_planar    = build_tensor(T, 'Saturated_Planar',    delta_lambdas, distance_lambdas) > 0.5;
sat_spherical = build_tensor(T, 'Saturated_Spherical', delta_lambdas, distance_lambdas) > 0.5;
sat_modified  = build_tensor(T, 'Saturated_Modified',  delta_lambdas, distance_lambdas) > 0.5;
sat_presonly  = build_tensor(T, 'Saturated_PresOnly',  delta_lambdas, distance_lambdas) > 0.5;
velocity_fd_error = build_tensor(T, 'Velocity_FD_Error', delta_lambdas, distance_lambdas);

% Output folder = same as CSV input
[results_folder, ~, ~] = fileparts(csv_path);
if isempty(results_folder)
    results_folder = '.';
end


%% STYLING — THESIS QUALITY %%

% Make sure LaTeX renders for all figure text
set(groot, 'defaulttextinterpreter', 'latex');
set(groot, 'defaultaxesticklabelinterpreter', 'latex');
set(groot, 'defaultlegendinterpreter', 'latex');

col_planar    = [0.250, 0.150, 0.040];
col_spherical = [0.550, 0.300, 0.080];
col_modified  = [0.850, 0.450, 0.100];
col_presonly  = [1.000, 0.700, 0.350];

fn = 'Times New Roman';
fs_ax = 13;
fs_lab = 15;
fs_leg = 12;
fs_title = 14;
lw_data = 1.8;
ms_data = 7;


%% FIGURE 1 — CROSSOVER MAP (MAIN BODY) %%
% Bracket difference V_3 Spherical minus V_4, displayed on a symmetric
% log colour scale so the crossover boundary is visible directly from
% the colour map (no contour line needed). Red = V_4 wins (smaller
% bracket), blue = V_3 Spherical wins.

bracket_3S_lambda = bracket_spherical / lambda;
bracket_V4_lambda = bracket_presonly  / lambda;
both_saturated = sat_spherical & sat_presonly;
bracket_diff = bracket_3S_lambda - bracket_V4_lambda;
bracket_diff_masked = bracket_diff;
bracket_diff_masked(both_saturated) = NaN;

% Symmetric log transform: sign(x) * log10(1 + |x|/linthresh)
linthresh = 0.05;
symlog = @(x) sign(x) .* log10(1 + abs(x) ./ linthresh);

display_data = symlog(bracket_diff_masked);
cmax_data = max(abs(bracket_diff_masked(:)), [], 'omitnan');
cmax_disp = symlog(cmax_data);

% Cell edges for axes (so ticks align with grid block boundaries)
x_edges = compute_cell_edges(distance_lambdas);
y_edges = compute_cell_edges(delta_lambdas);

figure('Color', 'w', 'Position', [100 100 760 620]);
data_padded = [display_data, nan(num_delta, 1); nan(1, num_dist + 1)];
pc = pcolor(x_edges, y_edges, data_padded);
pc.EdgeColor = 'none';
axis xy;

colormap(gca, redblue_colormap(256));
clim([-cmax_disp, cmax_disp]);

% Build colour bar with ticks at meaningful real-data values
cb_real_ticks = [-5, -1, -0.2, 0, 0.2, 1, 5];
cb_real_ticks = cb_real_ticks(abs(cb_real_ticks) <= cmax_data * 1.01);
cb = colorbar;
cb.Ticks = symlog(cb_real_ticks);
cb.TickLabels = arrayfun(@(v) strip_zeros(sprintf('%.2f', v)), ...
    cb_real_ticks, 'UniformOutput', false);
cb.Label.String = '$-3\,\mathrm{dB}$ range bracket advantage ($\lambda$)';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = fs_lab;
cb.TickLabelInterpreter = 'latex';

% Shrink axes + colourbar vertically to make room for endpoint labels
ax = gca;
ax_pos = ax.Position;
ax.Position = [ax_pos(1), ax_pos(2) + 0.04, ax_pos(3), ax_pos(4) - 0.08];
cb.Position = [cb.Position(1), cb.Position(2) + 0.04, ...
               cb.Position(3), cb.Position(4) - 0.08];

% Endpoint annotations — placed AFTER the position shift so coords are correct
cb_pos = cb.Position;
annotation('textbox', ...
    [cb_pos(1) - 0.02, cb_pos(2) + cb_pos(4) - 0.01, ...
     cb_pos(3) + 0.08, 0.035], ...
    'String', '$\mathbf{V}_4$ optimal', ...
    'Interpreter', 'latex', 'FontName', fn, 'FontSize', 11, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
    'EdgeColor', 'none', 'FitBoxToText', 'off');
annotation('textbox', ...
    [cb_pos(1) - 0.02, cb_pos(2) - 0.037, ...
     cb_pos(3) + 0.08, 0.035], ...
    'String', '$\mathbf{V}_3^{S}$ optimal', ...
    'Interpreter', 'latex', 'FontName', fn, 'FontSize', 11, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
    'EdgeColor', 'none', 'FitBoxToText', 'off');

hold on;

% Operational boundary lines
h_vlow = yline(velocity_lower_bound, 'k:',  'LineWidth', 1.5);
h_vhi  = yline(velocity_upper_bound, 'k-.', 'LineWidth', 1.5);

% Spatial aliasing band
alias_diagonal = 0.5 / sqrt(2);
xl = xlim;
hatch_spacing = (alias_threshold - alias_diagonal) / 2;
n_hatch = ceil((xl(2) - xl(1) + (alias_threshold - alias_diagonal)) / hatch_spacing);
for k = 0:n_hatch
    x_start = xl(1) + k * hatch_spacing - (alias_threshold - alias_diagonal);
    x_end   = x_start + (alias_threshold - alias_diagonal);
    xs = [max(x_start, xl(1)), min(x_end, xl(2))];
    ys = [alias_diagonal + (xs(1) - x_start), ...
          alias_diagonal + (xs(2) - x_start)];
    ys = max(min(ys, alias_threshold), alias_diagonal);
    plot(xs, ys, '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 0.8, ...
        'HandleVisibility', 'off');
end

yline(alias_diagonal,  'k-', 'LineWidth', 1.5);
yline(alias_threshold, 'k-', 'LineWidth', 1.5);
hold off;

xlabel('Source distance $r/\lambda$', 'FontSize', fs_lab);
ylabel('Array spacing $\delta/\lambda$', 'FontSize', fs_lab);

% Tick marks at every cell edge (visual grid delimiters), tick LABELS
% at every test value (where the data actually sits). Use minor ticks
% on edges, major ticks on test values. MATLAB doesn't quite do this
% directly, so we set major ticks at test values and rely on the cell
% boundaries being visually clear from the pcolor edges.

% Choose which test values to label — labelling every one is too dense
% on the x-axis where values cluster near zero.
x_label_vals = [0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0];
y_label_vals = [0.06  0.12  0.18  0.24  0.30  0.36  0.42];  % all 9 are fine on the y-axis

x_label_strs = arrayfun(@(v) strip_zeros(sprintf('%.2f', v)), ...
    x_label_vals, 'UniformOutput', false);
y_label_strs = arrayfun(@(v) strip_zeros(sprintf('%.2f', v)), ...
    y_label_vals, 'UniformOutput', false);

set(gca, 'FontName', fn, 'FontSize', fs_ax, 'LineWidth', 1.0, ...
    'Box', 'on', 'TickLabelInterpreter', 'latex', ...
    'XTick', x_label_vals, 'XTickLabel', x_label_strs, ...
    'YTick', y_label_vals, 'YTickLabel', y_label_strs, ...
    'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'TickDir', 'out', 'TickLength', [0.012 0.012]);

% Place minor ticks on the cell edges (the visual grid boundaries)
ax = gca;
ax.XAxis.MinorTickValues = x_edges;
ax.YAxis.MinorTickValues = y_edges;

xlim([x_edges(1), x_edges(end)]);
ylim([y_edges(1), y_edges(end)]);


% Legend below the map. The aliasing-band swatch is built as a
% custom hatched patch object — we use a HATCHFILL-style approach by
% creating a real (off-screen) hatched patch as the legend proxy.
% Trick: legend can render a pattern fill via the 'PolyShape' overlay
% only for filled patches, so we instead place a small inset axes
% directly inside the legend box AFTER it's drawn.
h_band_proxy = patch(NaN, NaN, [1 1 1], ...
    'EdgeColor', 'k', 'LineWidth', 1.0);

lgd = legend([h_vlow, h_vhi, h_band_proxy], ...
    {'Velocity SNR floor', ...
     'Linearity breakdown', ...
     'Aliasing-onset band'}, ...
    'Orientation', 'horizontal', ...
    'Location', 'southoutside', 'FontSize', fs_leg, ...
    'Box', 'on', 'Interpreter', 'latex');

% Overlay diagonal hatching onto the aliasing-band swatch by querying
% the legend's actual icon positions via its hidden Children. This is
% more robust than guessing column widths.
drawnow;
fig = gcf;
% Force MATLAB to commit the legend layout
lgd_pos = lgd.Position;

% Find the third entry's icon (the patch). Legend internally stores
% icon objects; we walk through them by inspecting the Children of
% the legend axes after drawing.
% Robust approach: get the patch's data extent in legend coordinates.
% MATLAB stores legend icons in lgd.EntryContainer or via undocumented
% APIs; the safest cross-version method is to find the icon by
% iterating peer plot objects with the same Tag.
%
% Practical workaround: place the overlay axes by computing the 3rd-
% entry x position from the legend's known internal layout. Horizontal
% legends in MATLAB pack entries left-to-right with a swatch + text
% for each. We measure the legend text extent to find the 3rd swatch.
icon_w = 0.05;   % swatch width as fraction of figure
icon_h = lgd_pos(4) * 0.57;
% The 3rd entry sits roughly 2/3 of the way across the legend; offset
% by ~5% of legend width as a left-margin estimate.
swatch_x = lgd_pos(1) + lgd_pos(3) * 0.67575;
swatch_y = lgd_pos(2) + (lgd_pos(4) - icon_h) / 2;

ax_swatch = axes('Parent', fig, ...
    'Position', [swatch_x, swatch_y, icon_w, icon_h]);
hold(ax_swatch, 'on');
patch(ax_swatch, [0 1 1 0], [0 0 1 1], 'w', ...
    'EdgeColor', 'k', 'LineWidth', 1.0);
n_hatch_sw = 2;
for k = -n_hatch_sw:n_hatch_sw
    xs = [k/n_hatch_sw, k/n_hatch_sw + 1];
    ys = [0, 1];
    xs_clip = max(min(xs, 1), 0);
    ys_clip = ys + (xs_clip - xs);
    plot(ax_swatch, xs_clip, ys_clip, '-', ...
        'Color', [0.3 0.3 0.3], 'LineWidth', 0.6);
end
xlim(ax_swatch, [0 1]); ylim(ax_swatch, [0 1]);
axis(ax_swatch, 'off');

save_thesis_figure(gcf, results_folder, 'crossover_map_thesis');


%% FIGURE 2 — SATURATION BOUNDARY (MAIN BODY) %%
% Operational range limit per method as a function of array spacing.

saturation_r = nan(num_delta, num_methods);
sat_tensors = {sat_planar, sat_spherical, sat_modified, sat_presonly};
for di = 1:num_delta
    for m = 1:num_methods
        sat_idx = find(sat_tensors{m}(di, :), 1, 'first');
        if ~isempty(sat_idx)
            saturation_r(di, m) = distance_lambdas(sat_idx);
        else
            saturation_r(di, m) = distance_lambdas(end) * 1.1;
        end
    end
end

figure('Color', 'w', 'Position', [100 100 620 420]);
plot(delta_lambdas, saturation_r(:, 1), '-o', ...
    'Color', col_planar, 'LineWidth', lw_data, 'MarkerSize', ms_data, ...
    'MarkerFaceColor', col_planar, 'DisplayName', 'Planar $\mathbf{V}_3$');
hold on;
plot(delta_lambdas, saturation_r(:, 2), '-d', ...
    'Color', col_spherical, 'LineWidth', lw_data, 'MarkerSize', ms_data, ...
    'MarkerFaceColor', col_spherical, 'DisplayName', 'Spherical $\mathbf{V}_3$');
plot(delta_lambdas, saturation_r(:, 3), '-s', ...
    'Color', col_modified, 'LineWidth', lw_data, 'MarkerSize', ms_data, ...
    'MarkerFaceColor', col_modified, 'DisplayName', 'Modified $\mathbf{V}_6$');
plot(delta_lambdas, saturation_r(:, 4), '-^', ...
    'Color', col_presonly, 'LineWidth', lw_data, 'MarkerSize', ms_data, ...
    'MarkerFaceColor', col_presonly, 'DisplayName', 'Pressure-only $\mathbf{V}_4$');

% Boundary lines
xline(velocity_lower_bound, 'k:', 'LineWidth', 1.3, ...
    'DisplayName', 'Velocity SNR floor');
xline(velocity_upper_bound, 'k-.', 'LineWidth', 1.3, ...
    'DisplayName', 'Linearity breakdown');
xline(alias_threshold, 'k--', 'LineWidth', 1.3, ...
    'DisplayName', 'Spatial Nyquist');
hold off;

xlabel('Array spacing $\delta/\lambda$', 'FontSize', fs_lab);
ylabel('First saturation distance $r/\lambda$', 'FontSize', fs_lab);
legend('Location', 'northwest', 'FontSize', fs_leg, 'Box', 'on', ...
    'NumColumns', 2);
set(gca, 'FontName', fn, 'FontSize', fs_ax, ...
    'XTick', delta_lambdas, 'LineWidth', 1.0, 'Box', 'on', ...
    'YGrid', 'on', 'XGrid', 'off', 'TickLabelInterpreter', 'latex');
save_thesis_figure(gcf, results_folder, 'saturation_boundary_thesis');


%% FIGURE 3 — VELOCITY FINITE-DIFFERENCE ERROR MAP (MAIN BODY) %%
% Validates the §3.4 upper boundary delta/lambda <= 0.33. Should be
% small in the operational range and grow above it.

figure('Color', 'w', 'Position', [100 100 700 480]);
% Cap colour scale at 1.0 — anything above means complete breakdown
% and using one limit keeps the colour map readable
display_max = 1.0;
display_data = velocity_fd_error;
display_data(display_data > display_max) = display_max;

imagesc(distance_lambdas, delta_lambdas, display_data);
axis xy;
colormap(gca, parula);
clim([0, display_max]);
cb = colorbar;
cb.Label.String = 'Relative finite-difference error';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = fs_lab;
cb.TickLabelInterpreter = 'latex';
hold on;

% Boundary lines in white for visibility on parula colour map
yline(velocity_lower_bound, 'w:', 'LineWidth', 1.8);
yline(velocity_upper_bound, 'w-.', 'LineWidth', 1.8);
yline(alias_threshold, 'w--', 'LineWidth', 1.8);

% Boundary labels
xl = xlim;
text(xl(2)*0.98, velocity_lower_bound, 'Velocity SNR floor', ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
    'FontName', fn, 'FontSize', 10, 'Color', 'w', 'Interpreter', 'latex');
text(xl(2)*0.98, velocity_upper_bound, 'Linearity breakdown', ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
    'FontName', fn, 'FontSize', 10, 'Color', 'w', 'Interpreter', 'latex');
text(xl(2)*0.98, alias_threshold, 'Spatial Nyquist', ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
    'FontName', fn, 'FontSize', 10, 'Color', 'w', 'Interpreter', 'latex');
hold off;

xlabel('Source distance $r/\lambda$', 'FontSize', fs_lab);
ylabel('Array spacing $\delta/\lambda$', 'FontSize', fs_lab);
set(gca, 'FontName', fn, 'FontSize', fs_ax, 'LineWidth', 1.0, ...
    'Box', 'on', 'TickLabelInterpreter', 'latex');
save_thesis_figure(gcf, results_folder, 'velocity_error_map_thesis');


%% FIGURE 4 — BRACKET GRID (APPENDIX) %%
% Family of bracket curves, one subplot per delta/lambda. Shows the
% raw data underlying the crossover map. Less synthesised than the
% main-body figures, suitable for appendix completeness.

n_subplot_rows = ceil(sqrt(num_delta));
n_subplot_cols = ceil(num_delta / n_subplot_rows);

figure('Color', 'w', 'Position', [100 100 1200 850]);

for di = 1:num_delta
    subplot(n_subplot_rows, n_subplot_cols, di);

    semilogy(distance_lambdas, bracket_planar(di, :)    / lambda, '-o', ...
        'Color', col_planar, 'LineWidth', 1.4, 'MarkerSize', 5, ...
        'MarkerFaceColor', col_planar);
    hold on;
    semilogy(distance_lambdas, bracket_spherical(di, :) / lambda, '-d', ...
        'Color', col_spherical, 'LineWidth', 1.4, 'MarkerSize', 5, ...
        'MarkerFaceColor', col_spherical);
    semilogy(distance_lambdas, bracket_modified(di, :)  / lambda, '-s', ...
        'Color', col_modified, 'LineWidth', 1.4, 'MarkerSize', 5, ...
        'MarkerFaceColor', col_modified);
    semilogy(distance_lambdas, bracket_presonly(di, :)  / lambda, '-^', ...
        'Color', col_presonly, 'LineWidth', 1.4, 'MarkerSize', 5, ...
        'MarkerFaceColor', col_presonly);
    hold off;

    title(sprintf('$\\delta/\\lambda = %.2f$', delta_lambdas(di)), ...
        'FontSize', 12);

    if di > num_delta - n_subplot_cols
        xlabel('$r/\lambda$', 'FontSize', 12);
    end
    if mod(di - 1, n_subplot_cols) == 0
        ylabel('$\Delta r/\lambda$', 'FontSize', 12);
    end
    if di == 1
        legend({'Planar $\mathbf{V}_3$', 'Spherical $\mathbf{V}_3$', ...
                'Modified $\mathbf{V}_6$', 'Pressure-only $\mathbf{V}_4$'}, ...
            'Location', 'southeast', 'FontSize', 9, 'Box', 'on');
    end

    set(gca, 'FontName', fn, 'FontSize', 10, ...
        'LineWidth', 0.8, 'Box', 'on', 'YGrid', 'on', ...
        'TickLabelInterpreter', 'latex');

    % Regime annotation
    dl = delta_lambdas(di);
    if dl < velocity_lower_bound
        regime_label = 'Below V-SNR floor';
        regime_colour = [0.7 0 0.7];
    elseif dl > alias_threshold
        regime_label = 'Aliased';
        regime_colour = 'r';
    elseif dl > velocity_upper_bound
        regime_label = 'V$_4$-only zone';
        regime_colour = [0.85 0.45 0.10];
    else
        regime_label = '';
    end
    if ~isempty(regime_label)
        text(0.5, 0.94, regime_label, 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'Color', regime_colour, ...
            'FontWeight', 'bold', 'FontSize', 10, 'Interpreter', 'latex');
    end
end

sgtitle('$-3\,\mathrm{dB}$ range bracket by $\delta/\lambda$', ...
    'FontSize', fs_title);
save_thesis_figure(gcf, results_folder, 'bracket_grid_thesis');


fprintf('\nAll thesis figures saved to: %s\n', results_folder);


%% LOCAL FUNCTIONS %%

function tensor = build_tensor(T, col_name, delta_lambdas, distance_lambdas)
    % Convert long-form CSV column into [num_delta x num_dist] tensor
    num_delta = length(delta_lambdas);
    num_dist = length(distance_lambdas);
    tensor = nan(num_delta, num_dist);
    col_data = T.(col_name);
    for di = 1:num_delta
        for ri = 1:num_dist
            mask = (T.Delta_Lambda == delta_lambdas(di)) & ...
                   (T.Distance_Lambda == distance_lambdas(ri));
            if any(mask)
                tensor(di, ri) = col_data(find(mask, 1));
            end
        end
    end
end


function save_thesis_figure(fig_handle, results_folder, base_name)
    exportgraphics(fig_handle, fullfile(results_folder, [base_name, '.png']), ...
        'ContentType', 'image', 'Resolution', 300);
    exportgraphics(fig_handle, fullfile(results_folder, [base_name, '.pdf']), ...
        'ContentType', 'vector');
    savefig(fig_handle, fullfile(results_folder, [base_name, '.fig']));
end


function cmap = redblue_colormap(N)
    % Diverging blue-white-orange map, neutral at zero. Anchor colours
    % chosen for thesis-quality contrast: deep blue at the negative
    % end, white at zero, burnt orange at the positive end.
    blue   = [0.10, 0.10, 1.00];
    white  = [1.00, 1.00, 1.00];
    orange = [0.9, 0.4, 0];

    if mod(N, 2) == 0
        half = N / 2;
        bot = [linspace(blue(1),  white(1), half)', ...
               linspace(blue(2),  white(2), half)', ...
               linspace(blue(3),  white(3), half)'];
        top = [linspace(white(1), orange(1), half)', ...
               linspace(white(2), orange(2), half)', ...
               linspace(white(3), orange(3), half)'];
    else
        half = (N - 1) / 2;
        bot = [linspace(blue(1),  white(1), half)', ...
               linspace(blue(2),  white(2), half)', ...
               linspace(blue(3),  white(3), half)'];
        top = [linspace(white(1), orange(1), half)', ...
               linspace(white(2), orange(2), half)', ...
               linspace(white(3), orange(3), half)'];
        bot = [bot; white];
    end
    cmap = [bot; top];
end


function edges = compute_cell_edges(test_values)
    % Returns the n+1 cell boundary positions for n test values.
    % Interior edges = midpoints between adjacent test values.
    % Outer edges = extrapolated by the same half-step.
    v = test_values(:)';
    mids = (v(1:end-1) + v(2:end)) / 2;
    left_edge  = v(1)   - (mids(1)   - v(1));
    right_edge = v(end) + (v(end)    - mids(end));
    edges = [left_edge, mids, right_edge];
end

function labels = format_edge_labels(test_values, edges)
    % For each edge position, find the nearest test value and label it.
    % Edges that don't correspond closely to a test value get an empty label
    % (this trims the two outer boundary edges which have no test value).
    n_edges = length(edges);
    labels  = cell(1, n_edges);
    for i = 1:n_edges
        [dist_to_nearest, idx] = min(abs(test_values - edges(i)));
        cell_spacing = median(diff(test_values));
        if dist_to_nearest < cell_spacing * 0.6
            v = test_values(idx);
            % Use 2 decimal places; strip trailing zeros for cleanliness
            s = sprintf('%.2f', v);
            s = regexprep(s, '\.?0+$', '');
            labels{i} = s;
        else
            labels{i} = '';  % outer boundary edges — no label
        end
    end
end

function s = strip_zeros(s)
    % Remove trailing decimal zeros: '0.50' -> '0.5', '1.00' -> '1'
    s = regexprep(s, '(\.\d*?)0+$', '$1');
    s = regexprep(s, '\.$', '');
end

function [tick_pos, tick_lab] = edge_ticks(test_values, edges)
    % Return tick positions at every cell edge, with labels only on
    % edges that coincide with a test value (within 5% of the LOCAL
    % cell spacing — handles non-uniform grids).
    tick_pos = edges;
    tick_lab = cell(1, length(edges));
    for i = 1:length(edges)
        [d, idx] = min(abs(test_values - edges(i)));
        % Use the local spacing around this test value
        if idx == 1
            local_spacing = test_values(2) - test_values(1);
        elseif idx == length(test_values)
            local_spacing = test_values(end) - test_values(end-1);
        else
            local_spacing = min(test_values(idx) - test_values(idx-1), ...
                                test_values(idx+1) - test_values(idx));
        end
        if d < local_spacing * 0.05
            tick_lab{i} = strip_zeros(sprintf('%.2f', test_values(idx)));
        else
            tick_lab{i} = '';
        end
    end
end