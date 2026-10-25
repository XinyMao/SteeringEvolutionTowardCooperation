function fig3_finite_revised_matlab()
%FIG3_FINITE_REVISED_MATLAB Recreate the finalized editable Figure 3.
%
% This function reads the frozen processed data and does not run simulations.
% It writes:
%   figures/Figure3_finite_revised_matlab.fig  (native editable MATLAB figure)
%   figures/Figure3_finite_revised_matlab.pdf  (editable vector PDF)
%   figures/Figure3_finite_revised_matlab.svg  (editable vector SVG, if supported)
%   figures/Figure3_finite_revised_matlab.png  (600-dpi raster preview)
%
% From the repository root, run with:
%   addpath(fullfile(pwd, 'scripts'));
%   fig3_finite_revised_matlab;

%% Configuration
rootDir = fileparts(fileparts(mfilename('fullpath')));
derivedDir = fullfile(rootDir, 'figures', 'derived');
revisedDir = fullfile(derivedDir, 'revised_fig3');
outputStem = fullfile(rootDir, 'figures', 'Figure3_finite_revised_matlab');

cfg.k = 6;
cfg.n = 4;
cfg.M = 20;
cfg.omega = 0.003;
cfg.c = 1.0;
cfg.rPanelA = 2.0;
cfg.pPanelA = [0.002, 0.010, 0.015, 0.020, 0.030];
cfg.pPanelB = [0.005, 0.010, 0.020, 0.030];
cfg.theoryPoints = 700;
cfg.pngDpi = 600;
cfg.tolerance = 1e-12;

panelAPath = fullfile(derivedDir, 'fig3a_drift_curve.csv');
panelBPath = fullfile(revisedDir, 'fig3b_local_threshold_revised.csv');
panelCEmpiricalPath = fullfile(revisedDir, 'fig3c_selected_high_x_estimator.csv');
panelCTheoryPath = fullfile(revisedDir, 'fig3c_theory_revised.csv');
requiredFiles = {panelAPath, panelBPath, panelCEmpiricalPath, panelCTheoryPath};
for i = 1:numel(requiredFiles)
    assert(isfile(requiredFiles{i}), 'Missing frozen input: %s', requiredFiles{i});
end

panelA = readtable(panelAPath);
panelB = readtable(panelBPath);
panelCEmpirical = readtable(panelCEmpiricalPath);
panelCTheory = readtable(panelCTheoryPath);

colorsA = [49 95 138; 111 159 193; 85 160 149; 226 140 69; 168 88 124] / 255;
colorsB = colorsA(1:4, :);
markers = {'o', 's', '^', 'd', 'v'};

%% Figure canvas
fig = figure('Color', 'w', 'Units', 'inches', ...
    'Position', [0.5 0.5 7.35 2.75], ...
    'Renderer', 'painters', ...
    'Name', 'Figure 3 finite steering', ...
    'NumberTitle', 'off');

% Use independent axes with explicit normalized positions.  This is
% intentional: tiledlayout constrains the axes as a layout group and makes
% manual frame dragging less convenient in MATLAB Plot Edit mode.
axA = axes('Parent', fig, 'Units', 'normalized', ...
    'Position', [0.070 0.190 0.255 0.685], 'Color', 'w');
axB = axes('Parent', fig, 'Units', 'normalized', ...
    'Position', [0.382 0.190 0.255 0.685], 'Color', 'w');
axC = axes('Parent', fig, 'Units', 'normalized', ...
    'Position', [0.694 0.190 0.255 0.685], 'Color', 'w');

%% Panel a: finite short-time drift
hold(axA, 'on');
for i = 1:numel(cfg.pPanelA)
    p = cfg.pPanelA(i);
    q = panelA(abs(panelA.p - p) < 1e-10, :);
    q = sortrows(q, 'x');
    assert(all(q.x >= xminFinite(p, cfg.k) - cfg.tolerance), ...
        'Panel-a data outside the corrected physical domain for p=%g.', p);

    % theory_drift is the frozen corrected finite-theory value used in the final plot.
    plot(axA, q.x, q.theory_drift, '-', 'Color', colorsA(i, :), ...
        'LineWidth', 1.35, 'HandleVisibility', 'off');
    neg = q.empirical_drift - q.CI_low;
    pos = q.CI_high - q.empirical_drift;
    errorbar(axA, q.x, q.empirical_drift, neg, pos, ...
        'LineStyle', 'none', 'Marker', markers{i}, 'MarkerSize', 4.2, ...
        'MarkerFaceColor', 'w', 'MarkerEdgeColor', colorsA(i, :), ...
        'Color', colorsA(i, :), 'LineWidth', 0.75, 'CapSize', 4, ...
        'DisplayName', sprintf('$p=%g$', p));
end
yline(axA, 0, '--', 'Color', [0.55 0.55 0.55], 'LineWidth', 0.8, ...
    'HandleVisibility', 'off');
xlim(axA, [0.39 0.975]);
xlabel(axA, 'Free-cooperator density, $x$', 'Interpreter', 'latex');
ylabel(axA, 'Short-time drift, $\langle\dot{x}\rangle$', 'Interpreter', 'latex');
legend(axA, 'Location', 'northeast', 'NumColumns', 2, 'Interpreter', 'latex', ...
    'FontSize', 7, 'Box', 'off');
formatAxes(axA);
panelLabel(axA, 'a');

%% Panel b: local zero-drift boundary
hold(axB, 'on');
pHandles = gobjects(numel(cfg.pPanelB), 1);
for i = 1:numel(cfg.pPanelB)
    p = cfg.pPanelB(i);
    x0 = xminFinite(p, cfg.k);
    xCurve = linspace(x0 + 1e-8, 1 - 1e-5, cfg.theoryPoints)';
    finiteCurve = rstarFinite(xCurve, p, cfg);
    sparseCurve = rstarSparse(xCurve, p, cfg);
    plot(axB, xCurve, finiteCurve, '-', 'Color', colorsB(i, :), ...
        'LineWidth', 1.35, 'HandleVisibility', 'off');
    plot(axB, xCurve, sparseCurve, '--', 'Color', colorsB(i, :), ...
        'LineWidth', 1.0, 'HandleVisibility', 'off');

    mask = abs(panelB.p - p) < 1e-10 & truthMask(panelB.bracketed);
    q = sortrows(panelB(mask, :), 'x');
    neg = q.r_empirical - q.CI_low;
    pos = q.CI_high - q.r_empirical;
    pHandles(i) = errorbar(axB, q.x, q.r_empirical, neg, pos, ...
        'LineStyle', 'none', 'Marker', markers{i}, 'MarkerSize', 4.4, ...
        'MarkerFaceColor', 'w', 'MarkerEdgeColor', colorsB(i, :), ...
        'Color', colorsB(i, :), 'LineWidth', 0.75, 'CapSize', 4, ...
        'DisplayName', sprintf('$p=%g$', p));

    weak = ~truthMask(q.envelope_eligible);
    if any(weak)
        plot(axB, q.x(weak), q.r_empirical(weak), 'x', ...
            'Color', colorsB(i, :), 'MarkerSize', 6, 'LineWidth', 1.0, ...
            'HandleVisibility', 'off');
    end
end
yline(axB, r0Finite(cfg.k, cfg.n), '--', 'Color', [0.55 0.55 0.55], ...
    'LineWidth', 0.8, 'HandleVisibility', 'off');
xlim(axB, [0.39 1.0]);
ylim(axB, [0.95 2.82]);
xlabel(axB, 'Free-cooperator density, $x$', 'Interpreter', 'latex');
ylabel(axB, 'Local threshold, $r^*(x,p)$', 'Interpreter', 'latex');
legend(axB, pHandles, 'Location', 'southeast', 'NumColumns', 2, ...
    'Interpreter', 'latex', 'FontSize', 6.7, 'Box', 'off');

% Editable manual style key, allowing the p legend to remain at lower right.
plot(axB, [0.415 0.475], [2.73 2.73], '-', 'Color', [0.15 0.15 0.15], ...
    'LineWidth', 1.35, 'HandleVisibility', 'off');
text(axB, 0.49, 2.73, 'finite theory', 'FontSize', 6.7, ...
    'VerticalAlignment', 'middle', 'Interpreter', 'none');
plot(axB, [0.415 0.475], [2.66 2.66], '--', 'Color', [0.15 0.15 0.15], ...
    'LineWidth', 1.0, 'HandleVisibility', 'off');
text(axB, 0.49, 2.66, 'sparse-limit baseline', 'FontSize', 6.7, ...
    'VerticalAlignment', 'middle', 'Interpreter', 'none');
formatAxes(axB);
panelLabel(axB, 'b');

%% Panel c: global analytical limit and high-x zero-drift estimates
hold(axC, 'on');
plot(axC, panelCTheory.p, panelCTheory.r_global_finite, 'k-', ...
    'LineWidth', 1.6, 'DisplayName', 'finite $r_g^*(p)$');
plot(axC, panelCTheory.p, panelCTheory.r_global_sparse, '--', ...
    'Color', [214 107 50] / 255, 'LineWidth', 1.35, ...
    'DisplayName', 'sparse-limit baseline');
neg = panelCEmpirical.estimate - panelCEmpirical.CI_low;
pos = panelCEmpirical.CI_high - panelCEmpirical.estimate;
errorbar(axC, panelCEmpirical.p, panelCEmpirical.estimate, neg, pos, ...
    'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', 5.0, ...
    'MarkerFaceColor', 'w', 'MarkerEdgeColor', [38 127 120] / 255, ...
    'Color', [38 127 120] / 255, 'LineWidth', 0.9, 'CapSize', 5, ...
    'DisplayName', 'high-$x$ zero-drift estimate');
xlim(axC, [0.001 0.035]);
ylim(axC, [0.8 2.68]);
xlabel(axC, 'Steering fraction, $p$', 'Interpreter', 'latex');
ylabel(axC, 'Multiplication factor, $r$', 'Interpreter', 'latex');
legend(axC, 'Location', 'southwest', 'Interpreter', 'latex', ...
    'FontSize', 6.7, 'Box', 'off');
formatAxes(axC);
panelLabel(axC, 'c');

% The title is a standalone annotation rather than a tiledlayout title, so it
% can be selected and dragged independently in the saved FIG file.
annotation(fig, 'textbox', [0.060 0.910 0.880 0.055], ...
    'String', 'Spatial correlations reshape cooperation boundaries under finite steering.', ...
    'Units', 'normalized', 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', 'FontName', 'Arial', 'FontSize', 11, ...
    'FontWeight', 'normal', 'Interpreter', 'none', 'LineStyle', 'none', ...
    'FitBoxToText', 'off', 'BackgroundColor', 'none', ...
    'Tag', 'Figure3EditableTitle');

%% Editable and publication outputs
drawnow;
savefig(fig, [outputStem '.fig'], 'compact');
exportgraphics(fig, [outputStem '.pdf'], 'ContentType', 'vector', ...
    'BackgroundColor', 'white');
exportgraphics(fig, [outputStem '.png'], 'Resolution', cfg.pngDpi, ...
    'BackgroundColor', 'white');
try
    exportgraphics(fig, [outputStem '.svg'], 'ContentType', 'vector', ...
        'BackgroundColor', 'white');
catch svgError
    warning('SVG export was skipped by this MATLAB version: %s', svgError.message);
end

fprintf('Wrote editable figure: %s.fig\n', outputStem);
fprintf('Wrote vector PDF:     %s.pdf\n', outputStem);
fprintf('Wrote PNG preview:    %s.png\n', outputStem);
if isfile([outputStem '.svg'])
    fprintf('Wrote vector SVG:     %s.svg\n', outputStem);
end
end

%% Local plotting and theory helpers
function formatAxes(ax)
set(ax, 'FontName', 'Arial', 'FontSize', 8, 'LineWidth', 0.8, ...
    'TickDir', 'out', 'Box', 'on', 'Layer', 'top');
end

function panelLabel(ax, label)
text(ax, -0.14, 1.04, label, 'Units', 'normalized', ...
    'FontName', 'Arial', 'FontSize', 13, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', ...
    'Clipping', 'off');
end

function mask = truthMask(values)
if islogical(values)
    mask = values;
elseif isnumeric(values)
    mask = values ~= 0;
else
    mask = strcmpi(string(values), "true") | string(values) == "1";
end
end

function value = r0Finite(k, n)
assert(k > 2 && n > 1);
value = k * n / (k + n - 1);
end

function value = xminFinite(p, k)
assert(k > 3);
assert(p > 0 && p < (k - 2) / (k - 1));
pc = (k - 3) / (2 * (k - 1));
if p <= pc
    a = 4 * p * (p + 1 / (k - 1));
    value = a / (1 + a);
else
    value = (k - 1) * p / (k - 2);
end
end

function value = rstarFinite(x, p, cfg)
assert(p > 0 && p < (cfg.k - 2) / (cfg.k - 1));
assert(all(x >= xminFinite(p, cfg.k) - cfg.tolerance));
assert(all(x < 1));

disc = x .* (x - 4 * p * (p + 1 / (cfg.k - 1)) .* (1 - x));
assert(all(disc >= -cfg.tolerance), ...
    'Negative discriminant outside floating-point tolerance.');
disc(disc < 0) = 0; % Only round-off-sized negative values reach this line.
delta = sqrt(disc);
eta0 = p - (0.5 + p) .* x + delta / 2;
assert(all(eta0 >= -cfg.tolerance));
eta0(eta0 < 0) = 0;
eta1 = (1 - x) .* ((cfg.k - 2) / (cfg.k - 1) .* x + 2 * eta0 - p);
eta2 = (cfg.n - 1) * (cfg.k - 1) / (cfg.k + cfg.n - 1) .* eta0 .* ...
    (1 - (1 - p) / p .* eta0);
denominator = eta1 + eta2;
assert(all(denominator > cfg.tolerance));
value = r0Finite(cfg.k, cfg.n) .* ...
    (eta1 - eta0 / (cfg.omega * cfg.M * cfg.c)) ./ denominator;
end

function value = rstarSparse(x, p, cfg)
assert(p > 0 && all(x > 0));
epsilon = cfg.omega * cfg.M * cfg.c;
assert(epsilon > 0);
value = r0Finite(cfg.k, cfg.n) .* (1 - p ./ (epsilon .* x));
end
