function fig = plot_sparse_final_figure2(cfg)
%PLOT_SPARSE_FINAL_FIGURE2 Final Nature Communications Figure 2.
%
%   plot_sparse_final_figure2
%   plot_sparse_final_figure2(cfg)
%
% The function reads the newest prepared sparse full-data pooled/cell tables,
% applies the PGG validity condition r=R*k*n/(k+n-1)>1 group by group, and
% exports the main figure, stand-alone panels, plot data, statistics, and
% the excluded-group table. It never modifies
% simulation code or source data and never runs simulations.

if nargin < 1 || isempty(cfg)
    cfg = struct();
end
cfg = apply_cfg_defaults(cfg);
if cfg.close_existing_figures, close all; end

scriptPath = mfilename("fullpath");
scriptDir = fileparts(scriptPath);
projectRoot = fileparts(fileparts(scriptDir));
resultsRoot = fullfile(projectRoot, "results");
mainOutputDir = fullfile(resultsRoot, "exp03_sparse_combined_pooled", "figures");
if ~isfolder(mainOutputDir)
    mkdir(mainOutputDir);
end

[dataDir, paths] = find_latest_sparse_data(resultsRoot);
fprintf("Final Figure 2 data directory: %s\n", dataDir);
fprintf("  cells: %s\n", paths.cells);
fprintf("  pooled bootstrap: %s\n", paths.bootstrap);
if strlength(paths.checklist) > 0
    fprintf("  checklist: %s\n", paths.checklist);
end

data = load_final_data(paths);
gAudit = audit_ghat_definition(projectRoot);
style = publication_style(cfg);
rng(20260802, "twister");

[fullGroups, validGroups, excludedGroups] = prepare_formal_groups(data.bootstrap);
fullCells = data.cells(data.cells.analysis == "full", :);
validCells = retain_valid_cells(fullCells, validGroups);

requiredR = [0.60, 0.80, 0.90];
loadedR = unique(fullGroups.R).';
if ~all(ismembertol(requiredR, loadedR, 1e-10))
    error("SparseFinalFig2:MissingNewR", ...
        "Latest pooled data do not contain all required R=0.60,0.80,0.90.");
end

panelA = prepare_equation_collapse(validCells, validGroups, paths.cells, ...
    cfg.panel_a_bins);
panelB = prepare_local_boundaries(validCells, paths.cells, cfg.panel_b_R, ...
    cfg.panel_b_ylim);
panelC = prepare_global_threshold(validGroups, validCells, paths.bootstrap, ...
    cfg.panel_c_only_well_identified);
panelD = prepare_coefficient_validation(panelC.displayData, paths.bootstrap);

excludedPath = fullfile(mainOutputDir, "Fig2_excluded_groups.csv");
writetable(excludedGroups, excludedPath);
weakPath = fullfile(mainOutputDir, "Fig2_weakly_identified_groups.csv");
writetable(panelC.weaklyIdentified, weakPath);

plotData = [panelA.plotRows; panelB.plotRows; panelC.plotRows; panelD.plotRows];
plotDataPath = fullfile(mainOutputDir, "Fig2_plot_data.csv");
writetable(plotData, plotDataPath);

statisticsPath = fullfile(mainOutputDir, "Fig2_statistics.txt");
write_statistics(statisticsPath, panelA, panelB, panelC, panelD, ...
    fullGroups, validGroups, excludedGroups, paths, gAudit);

[fig, mainAxes] = create_main_figure(panelA, panelB, panelC, panelD, style);
if cfg.output_variant == "refined"
    mainPdf = string(fullfile(mainOutputDir, "Fig2_sparse_refined.pdf"));
    mainPng = string(fullfile(mainOutputDir, "Fig2_sparse_refined.png"));
    mainFig = string(fullfile(mainOutputDir, "Fig2_sparse_refined.fig"));
    panelPaths = [
        string(fullfile(mainOutputDir, "Fig2a_equation_collapse.pdf"))
        string(fullfile(mainOutputDir, "Fig2b_local_threshold.pdf"))
        string(fullfile(mainOutputDir, "Fig2c_global_threshold.pdf"))
    ];
else
    mainPdf = string(fullfile(mainOutputDir, "Fig2_sparse_final.pdf"));
    mainPng = string(fullfile(mainOutputDir, "Fig2_sparse_final.png"));
    mainFig = string(fullfile(mainOutputDir, "Fig2_sparse_final.fig"));
    panelPaths = [
        string(fullfile(mainOutputDir, "Fig2_panel_a.pdf"))
        string(fullfile(mainOutputDir, "Fig2_panel_b.pdf"))
        string(fullfile(mainOutputDir, "Fig2_panel_c.pdf"))
    ];
end
if cfg.export_pdf
    exportgraphics(fig, mainPdf, "ContentType", "vector", "BackgroundColor", "white");
    for index = 1:2
        exportgraphics(mainAxes(index), panelPaths(index), ...
            "ContentType", "vector", "BackgroundColor", "white");
    end
    export_panel_c_split(panelC, style, panelPaths(3));
    export_panel_d_versions(panelD, style, mainOutputDir);
    if cfg.output_variant == "refined"
        export_refined_supplementary(panelB, fullGroups, validGroups, ...
            panelC, style, mainOutputDir);
    end
end
if cfg.export_legacy_supplementary
    balancedGroups = prepare_balanced_groups(data.bootstrap);
    export_supplementary_figures(validCells, fullCells, panelB, ...
        fullGroups, validGroups, balancedGroups, style, mainOutputDir, ...
        cfg.force_legacy_supplementary, cfg.legacy_supplementary_indices);
end
if cfg.export_png
    exportgraphics(fig, mainPng, "Resolution", cfg.png_dpi, ...
        "BackgroundColor", "white");
end
if cfg.export_fig
    savefig(fig, mainFig);
end
if ~cfg.keep_figure_open
    close(fig);
else
    fig.Visible = "on";
    figure(fig);
    drawnow;
    if cfg.enable_plot_edit
        plotedit(fig, "on");
    end
    if cfg.open_property_inspector
        inspect(fig);
    end
end

fprintf("\nFigure 2:\n");
if cfg.export_pdf, fprintf("  %s\n", mainPdf); end
if cfg.export_png, fprintf("  %s\n", mainPng); end
if cfg.export_fig, fprintf("  %s\n", mainFig); end
if cfg.export_pdf
    fprintf("Panels:\n");
    for index = 1:numel(panelPaths), fprintf("  %s\n", panelPaths(index)); end
end
fprintf("Panel a: Slope=%.6f, 95%% CI=[%.6f,%.6f], Weighted R2=%.6f\n", ...
    panelA.beta, panelA.betaCI, panelA.weightedR2);
fprintf("Panel c: Theory within CI=%d/%d, RMSE=%.6f\n", ...
    panelC.coverage, panelC.groups, panelC.rmse);
fprintf("Panel d theory coverage: alpha=%d/%d, b=%d/%d, s=%d/%d\n", ...
    panelD.alphaCoverage, panelD.groups, panelD.bCoverage, panelD.groups, ...
    panelD.sCoverage, panelD.groups);
fprintf("Panel b clipped CIs=%d; weakly identified points=%d\n", ...
    panelB.clippedCount, height(panelB.weaklyIdentified));
fprintf("Panel c PGG-valid=%d; well-identified=%d; weakly identified=%d\n", ...
    panelC.totalGroups, panelC.groups, height(panelC.weaklyIdentified));
fprintf("Colors: R .60=%s, .70=%s, .85=%s, .90=%s; epsilon .005=%s, .01=%s\n", ...
    cfg.panel_b_colors(1), cfg.panel_b_colors(2), ...
    cfg.panel_b_colors(3), cfg.panel_b_colors(4), ...
    cfg.palette(1), cfg.palette(4));
fprintf("Markers: Panel b o,s,d,v; Panels c/d ^,p,>\n");
end


function cfg = apply_cfg_defaults(cfg)
    defaults.font = 'Arial';
    defaults.export_pdf = true;
    defaults.export_png = true;
    defaults.export_fig = false;
    defaults.figure_visible = false;
    defaults.keep_figure_open = false;
    defaults.close_existing_figures = false;
    defaults.enable_plot_edit = false;
    defaults.open_property_inspector = false;
    defaults.panel_a_alpha = 0.15;
    defaults.panel_a_bins = 15;
    defaults.clip_panel_b = true;
    defaults.clip_panel_c = true;
    defaults.show_colorbar = true;
    defaults.palette = ["#3A5F8A", "#79A7C8", "#63A69F", ...
        "#E4985D", "#B8BDC6"];
    defaults.panel_b_R = [0.60, 0.70, 0.85, 0.90];
    defaults.panel_b_colors = ["#0072B2", "#D55E00", "#009E73", "#CC79A7"];
    defaults.panel_b_ylim = [-0.12, 0.45];
    defaults.panel_c_xlim = [-0.02, 0.55];
    defaults.panel_c_ylim = [-0.10, 0.90];
    defaults.figure_width_cm = 18.8;
    defaults.figure_height_cm = 14.0;
    defaults.png_dpi = 600;
    defaults.cell_marker_size = 7;
    defaults.marker_size = 5.2;
    defaults.theory_line_width = 2.0;
    defaults.errorbar_width = 0.85;
    defaults.violin_alpha = 0.35;
    defaults.panel_d_style = "raincloud";
    defaults.panel_c_only_well_identified = true;
    defaults.output_variant = "final";
    defaults.export_legacy_supplementary = false;
    defaults.force_legacy_supplementary = false;
    defaults.export_legacy_supplementary_fig = true;
    defaults.legacy_supplementary_indices = 1:8;
    defaults.axis_font_size = 8;
    defaults.label_font_size = 9;
    defaults.panel_label_size = 11;
    defaults.annotation_font_size = 8;
    fields = fieldnames(defaults);
    for index = 1:numel(fields)
        name = fields{index};
        if ~isfield(cfg, name) || isempty(cfg.(name))
            cfg.(name) = defaults.(name);
        end
    end
    validateattributes(cfg.panel_a_alpha, {'numeric'}, ...
        {'scalar','>=',0,'<=',1}, mfilename, 'cfg.panel_a_alpha');
    validateattributes(cfg.panel_a_bins, {'numeric'}, ...
        {'scalar','integer','>=',2}, mfilename, 'cfg.panel_a_bins');
    if numel(cfg.palette) ~= 5
        error("SparseFinalFig2:Palette", "cfg.palette must contain five colors.");
    end
    if numel(cfg.panel_b_colors) ~= 4
        error("SparseFinalFig2:PanelBColors", ...
            "cfg.panel_b_colors must contain four colors.");
    end
    validateattributes(cfg.legacy_supplementary_indices, {'numeric'}, ...
        {'vector','integer','>=',1,'<=',8}, mfilename, ...
        'cfg.legacy_supplementary_indices');
end


function [dataDir, paths] = find_latest_sparse_data(resultsRoot)
    candidates = dir(fullfile(resultsRoot, "**", ...
        "sparse_exploratory_36groups_bootstrap.csv"));
    if isempty(candidates)
        error("SparseFinalFig2:NoLatestPooled", ...
            "No prepared 36-group sparse pooled bootstrap table was found.");
    end
    [~, order] = sort([candidates.datenum], "descend");
    candidates = candidates(order);
    selected = [];
    for index = 1:numel(candidates)
        candidateDir = candidates(index).folder;
        cellPath = fullfile(candidateDir, ...
            "sparse_exploratory_36groups_cells.csv");
        inventoryPath = fullfile(candidateDir, ...
            "sparse_exploratory_36groups_inventory.csv");
        if isfile(cellPath) && isfile(inventoryPath)
            selected = index;
            break;
        end
    end
    if isempty(selected)
        error("SparseFinalFig2:IncompleteLatestPooled", ...
            "Bootstrap candidates were found, but none has cells and inventory tables.");
    end
    dataDir = candidates(selected).folder;
    paths.bootstrap = string(fullfile(dataDir, ...
        "sparse_exploratory_36groups_bootstrap.csv"));
    paths.cells = string(fullfile(dataDir, ...
        "sparse_exploratory_36groups_cells.csv"));
    paths.inventory = string(fullfile(dataDir, ...
        "sparse_exploratory_36groups_inventory.csv"));
    paths.summary = string(fullfile(dataDir, ...
        "sparse_exploratory_36groups_R_summary.csv"));
    checklistFiles = dir(fullfile(resultsRoot, "**", "*CHECKLIST*.md"));
    if isempty(checklistFiles)
        paths.checklist = "";
    else
        [~, latest] = max([checklistFiles.datenum]);
        paths.checklist = string(fullfile(checklistFiles(latest).folder, ...
            checklistFiles(latest).name));
    end
end


function data = load_final_data(paths)
    data.bootstrap = readtable(paths.bootstrap, "TextType", "string", ...
        "VariableNamingRule", "preserve");
    data.cells = readtable(paths.cells, "TextType", "string", ...
        "VariableNamingRule", "preserve");
    data.inventory = readtable(paths.inventory, "TextType", "string", ...
        "VariableNamingRule", "preserve");
    if isfile(paths.summary)
        data.summary = readtable(paths.summary, "TextType", "string", ...
            "VariableNamingRule", "preserve");
    else
        data.summary = table();
    end
    booleanColumns = ["threshold_covered", "alpha_covered", ...
        "b_covered", "s_covered"];
    for column = booleanColumns
        data.bootstrap.(column) = parse_boolean_column( ...
            data.bootstrap.(column), column);
    end
    validate_columns(data.bootstrap, ["analysis", "k", "n", "M", ...
        "R", "epsilon", "P_star_theory", "P_star_hat_x1", ...
        "ci_low", "ci_high", "intercept_over_epsilon", ...
        "intercept_ci_low", "intercept_ci_high", "coef_p", ...
        "coef_p_ci_low", "coef_p_ci_high", "state_slope_hat", ...
        "state_slope_ci_low", "state_slope_ci_high"]);
    validate_columns(data.cells, ["analysis", "k", "n", "M", "R", ...
        "epsilon", "x0_target", "p_over_p_star", "num_accepted", ...
        "mean_x_after_relax", "mean_p_actual", "mean_P_actual", ...
        "mean_G_over_epsilon", "sem_G_over_epsilon"]);
end


function audit = audit_ghat_definition(projectRoot)
    runnerPath = fullfile(projectRoot, "scripts", "run_sparse_drift.py");
    if ~isfile(runnerPath)
        error("SparseFinalFig2:MissingRunnerAudit", ...
            "Cannot audit G_hat because runner is missing: %s", runnerPath);
    end
    runnerText = fileread(runnerPath);
    compact = regexprep(runnerText, "\s+", "");
    expected = "G_hat=(params.k-1.0)/(params.k-2.0)*drift/(1.0-x_after_relax)";
    audit.correctStructureFactor = contains(compact, expected);
    audit.runnerPath = string(runnerPath);
    audit.formula = ["G_hat=(delta_x/tau_mcs)/", ...
        "[((k-2)/(k-1))*(1-x_after_relax)]"];
    if ~audit.correctStructureFactor
        error("SparseFinalFig2:GhatAuditFailed", ...
            "Stored G_hat definition could not be verified; raw reconstruction is required.");
    end
    fprintf("G_hat normalization audit passed: full (k-2)/(k-1) structure factor present.\n");
end


function validate_columns(T, required)
    names = string(T.Properties.VariableNames);
    missing = required(~ismember(required, names));
    if ~isempty(missing)
        fprintf("Available variables:\n");
        fprintf("  %s\n", names);
        error("SparseFinalFig2:MissingColumns", ...
            "Missing columns: %s", strjoin(missing, ", "));
    end
end


function values = parse_boolean_column(rawValues, columnName)
    if islogical(rawValues)
        values = rawValues;
    elseif isnumeric(rawValues)
        values = rawValues ~= 0;
    else
        textValues = lower(strtrim(string(rawValues)));
        known = ismember(textValues, ["true", "false", "1", "0"]);
        if ~all(known)
            error("SparseFinalFig2:BooleanColumn", ...
                "Unrecognized values in %s.", columnName);
        end
        values = ismember(textValues, ["true", "1"]);
    end
end


function [full, valid, excluded] = prepare_formal_groups(bootstrap)
    full = bootstrap(bootstrap.analysis == "full", :);
    full = sortrows(full, ["R", "epsilon", "k", "n", "M"]);
    if height(full) ~= 36
        error("SparseFinalFig2:FormalCount", ...
            "Expected 36 full-data formal groups; found %d.", height(full));
    end
    full.r0_computed = full.k .* full.n ./ (full.k + full.n - 1);
    full.r_computed = full.R .* full.r0_computed;
    full.isPGGValid = full.r_computed > 1 + 1e-12;
    valid = full(full.isPGGValid, :);
    excluded = full(~full.isPGGValid, ...
        ["k", "n", "M", "R", "epsilon", "r0_computed", "r_computed"]);
    excluded.exclusion_reason = repmat("r <= 1; outside main-text PGG domain", ...
        height(excluded), 1);
    excluded.source_file = repmat("sparse_exploratory_36groups_bootstrap.csv", ...
        height(excluded), 1);
end


function balanced = prepare_balanced_groups(bootstrap)
    balanced = bootstrap(bootstrap.analysis == "balanced300", :);
    balanced = sortrows(balanced, ["R", "epsilon", "k", "n", "M"]);
    balanced.r0_computed = balanced.k .* balanced.n ./ ...
        (balanced.k + balanced.n - 1);
    balanced.r_computed = balanced.R .* balanced.r0_computed;
    balanced.isPGGValid = balanced.r_computed > 1 + 1e-12;
end


function validCells = retain_valid_cells(cells, validGroups)
    key = group_key(cells.k, cells.n, cells.M, cells.R, cells.epsilon);
    validKey = group_key(validGroups.k, validGroups.n, validGroups.M, ...
        validGroups.R, validGroups.epsilon);
    validCells = cells(ismember(key, validKey), :);
    if height(validCells) ~= height(validGroups) * 56
        error("SparseFinalFig2:ValidCellCount", ...
            "Expected %d valid cells; found %d.", ...
            height(validGroups) * 56, height(validCells));
    end
end


function key = group_key(k, n, M, R, epsilon)
    key = compose("%d|%d|%d|%.10f|%.12f", k, n, M, R, epsilon);
end


function style = publication_style(cfg)
    style.cfg = cfg;
    style.font = string(cfg.font);
    style.axisFont = cfg.axis_font_size;
    style.labelFont = cfg.label_font_size;
    style.smallFont = 7.0;
    style.annotationFont = cfg.annotation_font_size;
    style.panelFont = cfg.panel_label_size;
    style.axisWidth = 0.9;
    style.theoryWidth = cfg.theory_line_width;
    style.errorWidth = cfg.errorbar_width;
    style.markerSize = cfg.marker_size;
    style.palette = zeros(5, 3);
    for index = 1:5
        style.palette(index, :) = hex2rgb(cfg.palette(index));
    end
    style.gameMarkers = ["^", "p", ">"];
    style.epsilonColors = style.palette([1, 4], :);
    style.rColors = zeros(4, 3);
    for index = 1:4
        style.rColors(index, :) = hex2rgb(cfg.panel_b_colors(index));
    end
    style.rMarkers = ["o", "s", "d", "v"];
    anchorX = linspace(0, 1, 4);
    style.rColormap = interp1(anchorX, style.palette(1:4, :), ...
        linspace(0, 1, 256), "pchip");
    style.rColormap = min(max(style.rColormap, 0), 1);
end


function rgb = hex2rgb(hex)
    value = char(erase(string(hex), "#"));
    if numel(value) ~= 6
        error("SparseFinalFig2:HexColor", "Invalid hexadecimal color: %s", hex);
    end
    rgb = [hex2dec(value(1:2)), hex2dec(value(3:4)), ...
        hex2dec(value(5:6))] / 255;
end


function panel = prepare_equation_collapse(cells, validGroups, sourcePath, numberOfBins)
    T = cells;
    T.Z = T.mean_P_actual - T.mean_x_after_relax .* (1 - T.R);
    T.Y = T.mean_G_over_epsilon;
    T.weight = inverse_variance_weights(T.sem_G_over_epsilon, T.num_accepted);
    finite = isfinite(T.Z) & isfinite(T.Y) & isfinite(T.weight) & T.weight > 0;
    if any(~finite)
        warning("SparseFinalFig2:PanelANonfinite", ...
            "Excluded %d nonfinite panel-a cells.", sum(~finite));
    end
    T = T(finite, :);
    fit = weighted_linear_fit(T.Z, T.Y, T.weight);
    fitted = fit.intercept + fit.slope .* T.Z;
    panel.weightedRMSE = sqrt(sum(T.weight .* (T.Y - fitted).^2) / sum(T.weight));

    bootstrapReps = 1000;
    [groupId, ~] = findgroups(T.k, T.n, T.M, T.R, T.epsilon);
    groupCount = max(groupId);
    bootCoef = nan(bootstrapReps, 2);
    for bootstrapIndex = 1:bootstrapReps
        sampled = zeros(height(T), 1);
        position = 1;
        for groupIndex = 1:groupCount
            source = find(groupId == groupIndex);
            draw = source(randi(numel(source), numel(source), 1));
            sampled(position:position + numel(draw) - 1) = draw;
            position = position + numel(draw);
        end
        B = T(sampled, :);
        bootFit = weighted_linear_fit(B.Z, B.Y, B.weight);
        bootCoef(bootstrapIndex, :) = [bootFit.intercept, bootFit.slope];
    end
    panel.alpha = fit.intercept;
    panel.beta = fit.slope;
    panel.alphaCI = prctile(bootCoef(:, 1), [2.5, 97.5]);
    panel.betaCI = prctile(bootCoef(:, 2), [2.5, 97.5]);
    panel.weightedR2 = fit.rSquared;
    panel.cells = height(T);
    panel.accepted = sum(T.num_accepted);
    panel.bootstrapReps = bootstrapReps;

    [~, order] = sort(T.Z);
    binId = zeros(height(T), 1);
    binId(order) = min(numberOfBins, ...
        ceil((1:height(T)).' .* numberOfBins ./ height(T)));
    bins = table();
    for binIndex = 1:numberOfBins
        rows = find(binId == binIndex);
        zMean = weighted_mean(T.Z(rows), T.weight(rows));
        yMean = weighted_mean(T.Y(rows), T.weight(rows));
        bootMean = nan(bootstrapReps, 1);
        for bootstrapIndex = 1:bootstrapReps
            sample = rows(randi(numel(rows), numel(rows), 1));
            bootMean(bootstrapIndex) = weighted_mean(T.Y(sample), T.weight(sample));
        end
        ci = prctile(bootMean, [2.5, 97.5]);
        newRow = table(binIndex, numel(rows), zMean, yMean, ci(1), ci(2), ...
            'VariableNames', ["bin", "num_cells", "Z", "estimate", ...
            "lower_CI", "upper_CI"]);
        bins = [bins; newRow]; %#ok<AGROW>
    end
    panel.data = T;
    panel.bins = bins;
    panel.numberOfBins = numberOfBins;
    panel.plotRows = empty_plot_table();
    for row = 1:height(T)
        panel.plotRows = append_plot_row(panel.plotRows, "a_cell", ...
            T.k(row), T.n(row), T.M(row), T.R(row), T.epsilon(row), ...
            T.mean_x_after_relax(row), T.mean_P_actual(row), T.Y(row), ...
            T.Y(row) - 1.96 * T.sem_G_over_epsilon(row), ...
            T.Y(row) + 1.96 * T.sem_G_over_epsilon(row), T.Z(row), ...
            "cell", sourcePath);
    end
    for row = 1:height(bins)
        panel.plotRows = append_plot_row(panel.plotRows, "a_bin", ...
            NaN, NaN, NaN, NaN, NaN, NaN, bins.Z(row), ...
            bins.estimate(row), bins.lower_CI(row), bins.upper_CI(row), ...
            bins.Z(row), compose("%d-bin IVW", numberOfBins), sourcePath);
    end
    alphaCovered = panel.alphaCI(1) <= 0 && panel.alphaCI(2) >= 0;
    betaCovered = panel.betaCI(1) <= 1 && panel.betaCI(2) >= 1;
    if ~alphaCovered || ~betaCovered
        warning("SparseFinalFig2:EquationCollapseInconsistent", ...
            "Equation-level collapse remains inconsistent with theory.");
    end
    fprintf("\nPanel a equation collapse: groups=%d, cells=%d, accepted=%d\n", ...
        height(validGroups), panel.cells, panel.accepted);
    fprintf("  alpha=%.6f [%.6f,%.6f]; beta=%.6f [%.6f,%.6f]\n", ...
        panel.alpha, panel.alphaCI, panel.beta, panel.betaCI);
    fprintf("  weighted R2=%.6f; RMSE=%.6f\n", ...
        panel.weightedR2, panel.weightedRMSE);
end


function panel = prepare_local_boundaries(cells, sourcePath, displayR, plottingYLimits)
    rValues = [0.60, 0.70, 0.80, 0.85, 0.90];
    mask = cells.k == 6 & cells.n == 3 & cells.M == 15 & ...
        abs(cells.epsilon - 0.01) < 1e-12 & ...
        ismembertol(cells.R, rValues, 1e-10);
    T = cells(mask, :);
    rows = table();
    bootstrapReps = 1000;
    for rIndex = 1:numel(rValues)
        R = rValues(rIndex);
        xValues = sort(unique(T.x0_target(abs(T.R - R) < 1e-10))).';
        for x0 = xValues
            group = sortrows(T(abs(T.R - R) < 1e-10 & ...
                abs(T.x0_target - x0) < 1e-10, :), "mean_P_actual");
            distinctP = unique(round(group.mean_P_actual, 12));
            if numel(distinctP) < 3
                warning("SparseFinalFig2:LocalDistinctP", ...
                    "Skipped R=%.2f,x0=%.2f: only %d actual P values.", ...
                    R, x0, numel(distinctP));
                continue;
            end
            weights = inverse_variance_weights( ...
                group.sem_G_over_epsilon, group.num_accepted);
            fit = weighted_linear_fit(group.mean_P_actual, ...
                group.mean_G_over_epsilon, weights);
            crossing = -fit.intercept / fit.slope;
            bootCrossing = nan(bootstrapReps, 1);
            bootSlope = nan(bootstrapReps, 1);
            for bootstrapIndex = 1:bootstrapReps
                bootY = group.mean_G_over_epsilon + ...
                    group.sem_G_over_epsilon .* randn(height(group), 1);
                bootFit = weighted_linear_fit(group.mean_P_actual, bootY, weights);
                bootSlope(bootstrapIndex) = bootFit.slope;
                if abs(bootFit.slope) > 1e-10
                    bootCrossing(bootstrapIndex) = -bootFit.intercept / bootFit.slope;
                end
            end
            validBoot = bootCrossing(isfinite(bootCrossing));
            ci = prctile(validBoot, [2.5, 97.5]);
            slopeCI = prctile(bootSlope(isfinite(bootSlope)), [2.5, 97.5]);
            sampledPMin = min(group.mean_P_actual);
            sampledPMax = max(group.mean_P_actual);
            slopeContainsZero = slopeCI(1) <= 0 && slopeCI(2) >= 0;
            crossingOutsideSampledP = crossing < sampledPMin || crossing > sampledPMax;
            finiteFraction = sum(isfinite(bootCrossing)) / bootstrapReps;
            weakIdentified = slopeContainsZero || crossingOutsideSampledP || ...
                finiteFraction < 0.80 || any(~isfinite(ci)) || numel(distinctP) < 3;
            actualX = weighted_mean(group.mean_x_after_relax, group.num_accepted);
            theory = actualX * (1 - R);
            covered = ci(1) <= theory && ci(2) >= theory;
            newRow = table(6, 3, 15, R, 0.01, x0, actualX, crossing, ...
                ci(1), ci(2), theory, fit.rSquared, numel(distinctP), ...
                sum(group.num_accepted), sum(isfinite(bootCrossing)), covered, ...
                slopeCI(1), slopeCI(2), sampledPMin, sampledPMax, ...
                finiteFraction, weakIdentified, ...
                'VariableNames', ["k", "n", "M", "R", "epsilon", ...
                "x0_target", "x", "estimate", "lower_CI", "upper_CI", ...
                "theory", "local_fit_R2", "num_P_cells", "num_accepted", ...
                "bootstrap_valid", "covered", "slope_ci_low", ...
                "slope_ci_high", "sampled_P_min", "sampled_P_max", ...
                "finite_bootstrap_fraction", "weakly_identified"]);
            rows = [rows; newRow]; %#ok<AGROW>
        end
    end
    panel.data = sortrows(rows, ["R", "x"]);
    panel.data.weak_reason = strings(height(panel.data), 1);
    for row = 1:height(panel.data)
        reasons = strings(0, 1);
        if panel.data.slope_ci_low(row) <= 0 && panel.data.slope_ci_high(row) >= 0
            reasons(end + 1) = "slope CI contains zero"; %#ok<AGROW>
        end
        if panel.data.estimate(row) < panel.data.sampled_P_min(row) || ...
                panel.data.estimate(row) > panel.data.sampled_P_max(row)
            reasons(end + 1) = "crossing outside sampled P"; %#ok<AGROW>
        end
        if panel.data.finite_bootstrap_fraction(row) < 0.80
            reasons(end + 1) = "finite bootstrap fraction below 0.80"; %#ok<AGROW>
        end
        if ~isfinite(panel.data.lower_CI(row)) || ~isfinite(panel.data.upper_CI(row))
            reasons(end + 1) = "nonfinite CI"; %#ok<AGROW>
        end
        if panel.data.num_P_cells(row) < 3
            reasons(end + 1) = "fewer than 3 actual-P cells"; %#ok<AGROW>
        end
        panel.data.weak_reason(row) = strjoin(reasons, "; ");
    end
    panel.rValues = rValues;
    panel.displayR = displayR;
    requested = ismembertol(panel.data.R, displayR, 1e-10);
    panel.weaklyIdentified = panel.data(requested & panel.data.weakly_identified, :);
    panel.displayData = panel.data(requested & ~panel.data.weakly_identified, :);
    panel.clippedCount = sum(panel.displayData.lower_CI < plottingYLimits(1) | ...
        panel.displayData.upper_CI > plottingYLimits(2));
    panel.bootstrapMethod = "cell-mean normal bootstrap; actual P and relaxed x";
    panel.plotRows = empty_plot_table();
    for row = 1:height(panel.displayData)
        panel.plotRows = append_plot_row(panel.plotRows, "b", ...
            panel.displayData.k(row), panel.displayData.n(row), panel.displayData.M(row), ...
            panel.displayData.R(row), panel.displayData.epsilon(row), panel.displayData.x(row), ...
            NaN, panel.displayData.estimate(row), panel.displayData.lower_CI(row), ...
            panel.displayData.upper_CI(row), panel.displayData.theory(row), ...
            "local WLS crossing", sourcePath);
    end
    fprintf("\nPanel b local boundaries:\n");
    for R = rValues
        group = panel.data(abs(panel.data.R - R) < 1e-10, :);
        errorValues = group.estimate - group.theory;
        fprintf("  R=%.2f: n=%d, median R2=%.4f, RMSE=%.5f, " + ...
            "MAE=%.5f, coverage=%d/%d\n", R, height(group), ...
            median(group.local_fit_R2), sqrt(mean(errorValues.^2)), ...
            mean(abs(errorValues)), sum(group.covered), height(group));
    end
end


function panel = prepare_global_threshold(validGroups, validCells, sourcePath, onlyWellIdentified)
    T = validGroups;
    T.theory = 1 - T.R;
    T.covered = T.ci_low <= T.theory & T.ci_high >= T.theory;
    T.sampled_P_min = nan(height(T), 1);
    T.sampled_P_max = nan(height(T), 1);
    for row = 1:height(T)
        mask = validCells.k == T.k(row) & validCells.n == T.n(row) & ...
            validCells.M == T.M(row) & abs(validCells.R - T.R(row)) < 1e-10 & ...
            abs(validCells.epsilon - T.epsilon(row)) < 1e-12;
        T.sampled_P_min(row) = min(validCells.mean_P_actual(mask));
        T.sampled_P_max(row) = max(validCells.mean_P_actual(mask));
    end
    T.slope_ci_contains_zero = T.coef_p_ci_low <= 0 & T.coef_p_ci_high >= 0;
    T.crossing_outside_sampled_P = T.P_star_hat_x1 < T.sampled_P_min | ...
        T.P_star_hat_x1 > T.sampled_P_max;
    T.finite_bootstrap_fraction = T.bootstrap_valid ./ T.bootstrap_reps;
    T.ci_nonfinite = ~isfinite(T.ci_low) | ~isfinite(T.ci_high);
    T.too_few_cells = T.num_cells < 3;
    T.weakly_identified = T.slope_ci_contains_zero | ...
        T.crossing_outside_sampled_P | T.finite_bootstrap_fraction < 0.80 | ...
        T.ci_nonfinite | T.too_few_cells;
    T.weak_reason = strings(height(T), 1);
    for row = 1:height(T)
        reasons = strings(0, 1);
        if T.slope_ci_contains_zero(row), reasons(end + 1) = "slope CI contains zero"; end %#ok<AGROW>
        if T.crossing_outside_sampled_P(row), reasons(end + 1) = "crossing outside sampled P"; end %#ok<AGROW>
        if T.finite_bootstrap_fraction(row) < 0.80, reasons(end + 1) = "finite bootstrap fraction below 0.80"; end %#ok<AGROW>
        if T.ci_nonfinite(row), reasons(end + 1) = "nonfinite CI"; end %#ok<AGROW>
        if T.too_few_cells(row), reasons(end + 1) = "fewer than 3 cells"; end %#ok<AGROW>
        T.weak_reason(row) = strjoin(reasons, "; ");
    end
    panel.data = T;
    panel.totalGroups = height(T);
    panel.weaklyIdentified = T(T.weakly_identified, :);
    if onlyWellIdentified
        panel.displayData = T(~T.weakly_identified, :);
    else
        panel.displayData = T;
    end
    panel.groups = height(panel.displayData);
    panel.coverage = sum(panel.displayData.covered);
    displayErrors = panel.displayData.P_star_hat_x1 - panel.displayData.theory;
    panel.rmse = sqrt(mean(displayErrors.^2));
    panel.mae = mean(abs(displayErrors));
    panel.uncovered = panel.displayData(~panel.displayData.covered, :);
    panel.plotRows = empty_plot_table();
    for row = 1:height(panel.displayData)
        D = panel.displayData;
        panel.plotRows = append_plot_row(panel.plotRows, "c", ...
            D.k(row), D.n(row), D.M(row), D.R(row), D.epsilon(row), ...
            NaN, D.theory(row), D.P_star_hat_x1(row), D.ci_low(row), ...
            D.ci_high(row), D.theory(row), "full pooled bootstrap", sourcePath);
    end
end


function panel = prepare_coefficient_validation(validGroups, sourcePath)
    T = validGroups;
    T.delta_alpha = T.intercept_over_epsilon;
    T.delta_alpha_low = T.intercept_ci_low;
    T.delta_alpha_high = T.intercept_ci_high;
    T.delta_b = T.coef_p - 1;
    T.delta_b_low = T.coef_p_ci_low - 1;
    T.delta_b_high = T.coef_p_ci_high - 1;
    T.delta_s = T.state_slope_hat - (1 - T.R);
    T.delta_s_low = T.state_slope_ci_low - (1 - T.R);
    T.delta_s_high = T.state_slope_ci_high - (1 - T.R);
    T.alphaCoveredFinal = T.delta_alpha_low <= 0 & T.delta_alpha_high >= 0;
    T.bCoveredFinal = T.delta_b_low <= 0 & T.delta_b_high >= 0;
    T.sCoveredFinal = T.delta_s_low <= 0 & T.delta_s_high >= 0;
    panel.data = T;
    panel.groups = height(T);
    panel.alphaCoverage = sum(T.alphaCoveredFinal);
    panel.bCoverage = sum(T.bCoveredFinal);
    panel.sCoverage = sum(T.sCoveredFinal);
    panel.summary = coefficient_summaries(T);
    panel.plotRows = empty_plot_table();
    categories = ["alpha", "b_minus_1", "s_minus_theory"];
    estimates = {T.delta_alpha, T.delta_b, T.delta_s};
    lows = {T.delta_alpha_low, T.delta_b_low, T.delta_s_low};
    highs = {T.delta_alpha_high, T.delta_b_high, T.delta_s_high};
    for categoryIndex = 1:3
        for row = 1:height(T)
            panel.plotRows = append_plot_row(panel.plotRows, "d", ...
                T.k(row), T.n(row), T.M(row), T.R(row), T.epsilon(row), ...
                NaN, categoryIndex, estimates{categoryIndex}(row), ...
                lows{categoryIndex}(row), highs{categoryIndex}(row), 0, ...
                categories(categoryIndex), sourcePath);
        end
    end
end


function summary = coefficient_summaries(T)
    names = ["alpha", "b_minus_1", "s_minus_theory"];
    estimates = {T.delta_alpha, T.delta_b, T.delta_s};
    lows = {T.delta_alpha_low, T.delta_b_low, T.delta_s_low};
    highs = {T.delta_alpha_high, T.delta_b_high, T.delta_s_high};
    summary = table();
    for index = 1:3
        standardError = (highs{index} - lows{index}) ./ (2 * 1.96);
        weights = 1 ./ standardError.^2;
        finite = isfinite(estimates{index}) & isfinite(weights) & weights > 0;
        estimate = weighted_mean(estimates{index}(finite), weights(finite));
        se = sqrt(1 / sum(weights(finite)));
        row = table(names(index), estimate, estimate - 1.96 * se, ...
            estimate + 1.96 * se, sum(finite), ...
            'VariableNames', ["coefficient", "estimate", "lower_CI", ...
            "upper_CI", "num_groups"]);
        summary = [summary; row]; %#ok<AGROW>
    end
end


function weights = inverse_variance_weights(semValues, counts)
    weights = zeros(size(semValues));
    valid = isfinite(semValues) & semValues > 0;
    weights(valid) = 1 ./ semValues(valid).^2;
    fallback = ~valid & isfinite(counts) & counts > 0;
    weights(fallback) = counts(fallback);
end


function value = weighted_mean(values, weights)
    finite = isfinite(values) & isfinite(weights) & weights > 0;
    value = sum(values(finite) .* weights(finite)) / sum(weights(finite));
end


function fit = weighted_linear_fit(x, y, weights)
    finite = isfinite(x) & isfinite(y) & isfinite(weights) & weights > 0;
    x = x(finite);
    y = y(finite);
    weights = weights(finite);
    if numel(x) < 3
        error("SparseFinalFig2:TooFewFitPoints", ...
            "Weighted fit requires at least three points.");
    end
    X = [ones(numel(x), 1), x];
    matrix = X.' * (weights .* X);
    if rcond(matrix) < 1e-12
        error("SparseFinalFig2:IllConditionedFit", "Weighted fit is ill-conditioned.");
    end
    coefficient = matrix \ (X.' * (weights .* y));
    fitted = X * coefficient;
    yMean = weighted_mean(y, weights);
    residualSS = sum(weights .* (y - fitted).^2);
    totalSS = sum(weights .* (y - yMean).^2);
    fit.intercept = coefficient(1);
    fit.slope = coefficient(2);
    fit.rSquared = 1 - residualSS / totalSS;
end


function [fig, axesMain] = create_main_figure(panelA, panelB, panelC, panelD, style)
    % Keep the window hidden while exports run. Otherwise MATLAB is busy but
    % the visible Format-tab controls appear unresponsive to the user.
    visibility = "off";
    fig = figure("Color", "white", "Units", "centimeters", ...
        "Position", [1, 1, style.cfg.figure_width_cm, ...
        style.cfg.figure_height_cm], "PaperPositionMode", "auto", ...
        "Visible", visibility, "HandleVisibility", "on");
    % Explicit axes positions are intentional: unlike tiledlayout-managed
    % axes, these frames remain independently editable in Property Inspector.
    axA = axes(fig, "Units", "normalized", ...
        "Position", [0.070 0.575 0.335 0.370], "Tag", "panel_a");
    draw_panel_a(axA, panelA, style, true);
    axB = axes(fig, "Units", "normalized", ...
        "Position", [0.620 0.575 0.340 0.370], "Tag", "panel_b");
    draw_panel_b(axB, panelB, style, true);
    axC1 = axes(fig, "Units", "normalized", ...
        "Position", [0.070 0.315 0.365 0.165], "Tag", "panel_c_eps0005");
    draw_panel_c(axC1, panelC, style, true, false, 0.005);
    axC2 = axes(fig, "Units", "normalized", ...
        "Position", [0.070 0.075 0.365 0.165], "Tag", "panel_c_eps001");
    draw_panel_c(axC2, panelC, style, false, false, 0.01);
    axD1 = axes(fig, "Units", "normalized", ...
        "Position", [0.600 0.315 0.360 0.165], "Tag", "panel_d_eps0005");
    draw_panel_d(axD1, panelD, style, true, 0.005);
    axD2 = axes(fig, "Units", "normalized", ...
        "Position", [0.600 0.075 0.360 0.165], "Tag", "panel_d_eps001");
    draw_panel_d(axD2, panelD, style, false, 0.01);
    axesMain = [axA, axB, axC1, axC2, axD1, axD2];
    for ax = axesMain
        ax.PositionConstraint = "innerposition";
    end
end


function draw_panel_a(ax, panel, style, showLabel)
    hold(ax, "on");
    box(ax, "on");
    values = [panel.data.Z; panel.data.Y; panel.bins.lower_CI; panel.bins.upper_CI];
    limits = data_limits(values, 0.05);
    cells = scatter(ax, panel.data.Z, panel.data.Y, ...
        style.cfg.cell_marker_size, panel.data.R, "filled", ...
        "MarkerEdgeColor", "none", ...
        "MarkerFaceAlpha", style.cfg.panel_a_alpha);
    colormap(ax, style.rColormap);
    clim(ax, [0.5 0.9]);
    binned = errorbar(ax, panel.bins.Z, panel.bins.estimate, ...
        panel.bins.estimate - panel.bins.lower_CI, ...
        panel.bins.upper_CI - panel.bins.estimate, ...
        "LineStyle", "none", "Marker", "o", "MarkerSize", 5.8, ...
        "MarkerFaceColor", [0.08 0.08 0.08], ...
        "MarkerEdgeColor", "white", "Color", [0.08 0.08 0.08], ...
        "LineWidth", style.errorWidth, "CapSize", 2.5);
    theory = plot(ax, limits, limits, "k--", "LineWidth", style.theoryWidth);
    xlabel(ax, {"Predicted normalized drift,", "$Z=P-x(1-R)$"}, ...
        "Interpreter", "latex");
    ylabel(ax, "$\widehat{G}/\epsilon$", "Interpreter", "latex");
    legend(ax, [cells, binned, theory], ...
        {"Cell means", "Binned means (95% CI)", "Theory"}, ...
        "Location", "northwest", "Box", "off", "FontSize", style.smallFont);
    if style.cfg.show_colorbar
        colorScale = colorbar(ax);
        colorScale.Ticks = [0.5 0.7 0.9];
        colorScale.TickLabels = {"0.5", "0.7", "0.9"};
        colorScale.Label.String = "R";
        colorScale.FontName = style.font;
        colorScale.FontSize = style.axisFont;
    end
    text(ax, 0.97, 0.05, ...
        "Slope = " + compose("%.3f", panel.beta) + newline + ...
        "95% CI [" + compose("%.3f", panel.betaCI(1)) + ", " + ...
        compose("%.3f", panel.betaCI(2)) + "]" + newline + ...
        "Weighted R^2 = " + compose("%.2f", panel.weightedR2), ...
        "Units", "normalized", "HorizontalAlignment", "right", ...
        "VerticalAlignment", "bottom", ...
        "FontSize", style.annotationFont);
    xlim(ax, limits);
    ylim(ax, limits);
    apply_style(ax, style);
    if showLabel, add_panel_label(ax, "a", style); end
end


function draw_panel_b(ax, panel, style, showLabel)
    hold(ax, "on");
    box(ax, "on");
    handles = gobjects(numel(panel.displayR), 1);
    yLimits = style.cfg.panel_b_ylim;
    for index = 1:numel(panel.displayR)
        R = panel.displayR(index);
        group = panel.displayData(abs(panel.displayData.R - R) < 1e-10, :);
        theoryX = linspace(0.2, 0.8, 100);
        theoryColor = 0.18 + 0.82 * style.rColors(index, :);
        plot(ax, theoryX, theoryX .* (1 - R), "-", ...
            "Color", theoryColor, ...
            "LineWidth", style.theoryWidth, "HandleVisibility", "off");
        if isempty(group)
            handles(index) = plot(ax, NaN, NaN, style.rMarkers(index), ...
                "MarkerFaceColor", style.rColors(index, :), ...
                "MarkerEdgeColor", [0.12 0.12 0.12], "LineStyle", "none");
        else
            handles(index) = draw_clipped_errorbars(ax, group.x, group.estimate, ...
                group.lower_CI, group.upper_CI, style.rColors(index, :), ...
                style.rMarkers(index), yLimits, style.cfg.clip_panel_b, style);
        end
    end
    xlabel(ax, "Cooperator fraction, $x$", "Interpreter", "latex");
    ylabel(ax, "$\widehat{P}^{*}(x)$", "Interpreter", "latex");
    legend(ax, handles, cellstr(compose("R=%.2f", panel.displayR)), ...
        "Location", "northwest", "Box", "off", ...
        "FontSize", style.smallFont);
    xlim(ax, [0.18 0.82]);
    ylim(ax, yLimits);
    apply_style(ax, style);
    if showLabel, add_panel_label(ax, "b", style); end
end


function draw_panel_c(ax, panel, style, showLabel, includeInvalid, epsilonFilter)
    if nargin < 6, epsilonFilter = NaN; end
    hold(ax, "on");
    box(ax, "on");
    T = panel.displayData;
    if isfinite(epsilonFilter)
        T = T(abs(T.epsilon - epsilonFilter) < 1e-12, :);
    end
    xLimits = style.cfg.panel_c_xlim;
    yLimits = style.cfg.panel_c_ylim;
    theoryLimits = [max(xLimits(1), yLimits(1)), min(xLimits(2), yLimits(2))];
    plot(ax, theoryLimits, theoryLimits, "k--", "LineWidth", style.theoryWidth, ...
        "HandleVisibility", "off");
    games = sortrows(unique(T(:, ["k", "n", "M"]), "rows"), ...
        ["k", "n", "M"]);
    epsilonValues = sort(unique(T.epsilon));
    handles = gobjects(3, 1);
    for epsilonIndex = 1:numel(epsilonValues)
        colorIndex = 1 + (epsilonValues(epsilonIndex) > 0.0075);
        for gameIndex = 1:height(games)
            mask = abs(T.epsilon - epsilonValues(epsilonIndex)) < 1e-12 & ...
                T.k == games.k(gameIndex) & T.n == games.n(gameIndex) & ...
                T.M == games.M(gameIndex);
            group = T(mask, :);
            handle = draw_clipped_errorbars(ax, group.theory, ...
                group.P_star_hat_x1, group.ci_low, group.ci_high, ...
                style.epsilonColors(colorIndex, :), ...
                style.gameMarkers(gameIndex), yLimits, ...
                style.cfg.clip_panel_c, style);
            handle.HandleVisibility = "off";
        end
    end
    for index = 1:3
        handles(index) = plot(ax, NaN, NaN, style.gameMarkers(index), ...
            "MarkerSize", style.markerSize, "MarkerFaceColor", [0.65 0.65 0.65], ...
            "MarkerEdgeColor", [0.12 0.12 0.12], "LineStyle", "none");
    end
    if showLabel || ~isfinite(epsilonFilter)
        legend(ax, handles, {"(4,3,6)", "(4,4,4)", "(6,3,15)"}, ...
            "Location", "northwest", "NumColumns", 3, "Box", "off", ...
            "FontSize", style.smallFont);
    end
    if ~isfinite(epsilonFilter) || epsilonFilter > 0.0075
        xlabel(ax, "Theoretical global threshold, $1-R$", "Interpreter", "latex");
    else
        ax.XTickLabel = [];
    end
    ylabel(ax, "$\widehat{P}^{*}(1)$", "Interpreter", "latex");
    if ~includeInvalid
        localCoverage = sum(T.covered);
        localRMSE = sqrt(mean((T.P_star_hat_x1 - T.theory).^2));
        text(ax, 0.97, 0.06, sprintf("epsilon=%.3f; CI %d/%d; RMSE %.3f", ...
            epsilonValues(1), localCoverage, height(T), localRMSE), ...
            "Units", "normalized", "HorizontalAlignment", "right", ...
            "VerticalAlignment", "bottom", ...
            "FontSize", style.smallFont);
    end
    xlim(ax, xLimits);
    ylim(ax, yLimits);
    apply_style(ax, style);
    if showLabel, add_panel_label(ax, "c", style); end
end


function draw_panel_d(ax, panel, style, showLabel, epsilonFilter)
    if nargin < 5, epsilonFilter = NaN; end
    hold(ax, "on");
    box(ax, "on");
    yline(ax, 0, "k--", "LineWidth", 1.5, "HandleVisibility", "off");
    T = panel.data;
    if isfinite(epsilonFilter)
        T = T(abs(T.epsilon - epsilonFilter) < 1e-12, :);
    end
    games = sortrows(unique(T(:, ["k", "n", "M"]), "rows"), ...
        ["k", "n", "M"]);
    epsilonValues = sort(unique(T.epsilon));
    localSummary = coefficient_summaries(T);
    estimate = {T.delta_alpha, T.delta_b, T.delta_s};
    for category = 1:3
        if style.cfg.panel_d_style == "raincloud"
            for epsilonIndex = 1:numel(epsilonValues)
                colorIndex = 1 + (epsilonValues(epsilonIndex) > 0.0075);
                values = estimate{category}(abs(T.epsilon - ...
                    epsilonValues(epsilonIndex)) < 1e-12);
                draw_half_violin_right(ax, category + 0.08 + ...
                    0.025 * (epsilonIndex - 1), values, ...
                    style.epsilonColors(colorIndex, :), ...
                    style.cfg.violin_alpha, 0.22);
            end
        end
        offsets = deterministic_beeswarm(estimate{category}, 0.11);
        for epsilonIndex = 1:numel(epsilonValues)
            colorIndex = 1 + (epsilonValues(epsilonIndex) > 0.0075);
            for gameIndex = 1:3
                mask = abs(T.epsilon - epsilonValues(epsilonIndex)) < 1e-12 & ...
                    T.k == games.k(gameIndex) & T.n == games.n(gameIndex) & ...
                    T.M == games.M(gameIndex);
                rows = find(mask);
                scatter(ax, category - 0.10 + offsets(rows), ...
                    estimate{category}(rows), 18, ...
                    "Marker", style.gameMarkers(gameIndex), ...
                    "MarkerFaceColor", style.epsilonColors(colorIndex, :), ...
                    "MarkerEdgeColor", [0.15 0.15 0.15], ...
                    "MarkerFaceAlpha", 0.74, "MarkerEdgeAlpha", 0.80, ...
                    "HandleVisibility", "off");
            end
        end
        summary = localSummary(category, :);
        errorbar(ax, category + 0.02, summary.estimate, ...
            summary.estimate - summary.lower_CI, ...
            summary.upper_CI - summary.estimate, ...
            "LineStyle", "none", "Marker", "d", "MarkerSize", 7, ...
            "MarkerFaceColor", "black", "MarkerEdgeColor", "black", ...
            "Color", "black", "LineWidth", 1.6, "CapSize", 4, ...
            "HandleVisibility", "off");
    end
    ax.XTick = 1:3;
    ax.XTickLabel = {"$\widehat{\alpha}$", "$\widehat{b}-1$", ...
        "$\widehat{s}-(1-R)$"};
    ax.TickLabelInterpreter = "latex";
    if isfinite(epsilonFilter) && epsilonFilter < 0.0075
        ax.XTickLabel = [];
    end
    ylabel(ax, "Estimate - theory");
    xlim(ax, [0.40 4.00]);
    allT = panel.data;
    values = [allT.delta_alpha; allT.delta_b; allT.delta_s; ...
        panel.summary.lower_CI; panel.summary.upper_CI];
    ylim(ax, data_limits(values, 0.12));
    alphaCoverage = sum(T.delta_alpha_low <= 0 & T.delta_alpha_high >= 0);
    bCoverage = sum(T.delta_b_low <= 0 & T.delta_b_high >= 0);
    sCoverage = sum(T.delta_s_low <= 0 & T.delta_s_high >= 0);
    text(ax, 0.97, 0.95, sprintf("epsilon=%.3f; CI: alpha %d/%d, b %d/%d, s %d/%d", ...
        epsilonValues(1), alphaCoverage, height(T), bCoverage, height(T), ...
        sCoverage, height(T)), "Units", "normalized", ...
        "HorizontalAlignment", "right", "VerticalAlignment", "top", ...
        "FontSize", style.annotationFont);
    apply_style(ax, style);
    if showLabel, add_panel_label(ax, "d", style); end
end


function draw_half_violin_right(ax, category, values, color, alphaValue, maximumWidth)
    values = values(isfinite(values));
    if numel(values) < 2, return; end
    spread = std(values);
    if ~isfinite(spread) || spread <= 0
        spread = max(abs(values(1)), 1) * 0.05;
    end
    bandwidth = max(1.06 * spread * numel(values)^(-1/5), ...
        eps(max(abs(values)) + 1));
    yGrid = linspace(min(values) - 0.35 * spread, ...
        max(values) + 0.35 * spread, 160);
    density = zeros(size(yGrid));
    for index = 1:numel(values)
        density = density + exp(-0.5 * ((yGrid - values(index)) / bandwidth).^2);
    end
    if max(density) > 0, density = maximumWidth * density / max(density); end
    patch(ax, [category + density, category * ones(size(density))], ...
        [yGrid, fliplr(yGrid)], color, "FaceAlpha", alphaValue, ...
        "EdgeColor", "none", "HandleVisibility", "off");
end


function offsets = deterministic_beeswarm(values, maximumWidth)
    offsets = zeros(size(values));
    finiteRows = find(isfinite(values));
    if numel(finiteRows) < 2, return; end
    finiteValues = values(finiteRows);
    edges = linspace(min(finiteValues), max(finiteValues), 11);
    if edges(1) == edges(end), return; end
    bins = discretize(finiteValues, edges);
    bins(isnan(bins)) = 10;
    for bin = 1:10
        rows = finiteRows(bins == bin);
        [~, order] = sort(values(rows));
        rows = rows(order);
        if isscalar(rows)
            local = 0;
        else
            local = linspace(-maximumWidth, maximumWidth, numel(rows));
        end
        offsets(rows) = local;
    end
end


function handle = draw_clipped_errorbars(ax, x, estimate, lowerCI, upperCI, ...
        color, marker, yLimits, clipIntervals, style)
    if clipIntervals
        plottedLow = max(lowerCI, yLimits(1));
        plottedHigh = min(upperCI, yLimits(2));
        plottedEstimate = min(max(estimate, yLimits(1)), yLimits(2));
        plottedLow = min(plottedLow, plottedEstimate);
        plottedHigh = max(plottedHigh, plottedEstimate);
    else
        plottedLow = lowerCI;
        plottedHigh = upperCI;
        plottedEstimate = estimate;
    end
    handle = errorbar(ax, x, plottedEstimate, plottedEstimate - plottedLow, ...
        plottedHigh - plottedEstimate, "LineStyle", "none", "Marker", marker, ...
        "MarkerSize", style.markerSize, "MarkerFaceColor", color, ...
        "MarkerEdgeColor", [0.12 0.12 0.12], "Color", color, ...
        "LineWidth", style.errorWidth, "CapSize", 2.2);
    if clipIntervals
        span = diff(yLimits);
        lowerClipped = lowerCI < yLimits(1) | estimate < yLimits(1);
        upperClipped = upperCI > yLimits(2) | estimate > yLimits(2);
        plot(ax, x(lowerClipped), repmat(yLimits(1) + 0.012 * span, ...
            sum(lowerClipped), 1), "v", "Color", color, ...
            "MarkerFaceColor", color, "MarkerSize", 3.3, ...
            "LineStyle", "none", "HandleVisibility", "off");
        plot(ax, x(upperClipped), repmat(yLimits(2) - 0.012 * span, ...
            sum(upperClipped), 1), "^", "Color", color, ...
            "MarkerFaceColor", color, "MarkerSize", 3.3, ...
            "LineStyle", "none", "HandleVisibility", "off");
    end
end


function export_panel_d_versions(panel, style, outDir)
    names = ["Fig2d_coefficients_raincloud.pdf", ...
        "Fig2d_coefficients_beeswarm.pdf"];
    variants = ["raincloud", "beeswarm"];
    for index = 1:2
        localStyle = style;
        localStyle.cfg.panel_d_style = variants(index);
        fig = standalone_figure([9.2, 9.0]);
        axTop = axes(fig, "Units", "normalized", ...
            "Position", [0.14 0.56 0.80 0.36]);
        draw_panel_d(axTop, panel, localStyle, true, 0.005);
        axBottom = axes(fig, "Units", "normalized", ...
            "Position", [0.14 0.10 0.80 0.36]);
        draw_panel_d(axBottom, panel, localStyle, false, 0.01);
        exportgraphics(fig, fullfile(outDir, names(index)), ...
            "ContentType", "vector", "BackgroundColor", "white");
        close(fig);
    end
end


function export_panel_c_split(panel, style, path)
    fig = standalone_figure([9.2, 9.0]);
    axTop = axes(fig, "Units", "normalized", ...
        "Position", [0.14 0.56 0.80 0.36]);
    draw_panel_c(axTop, panel, style, true, false, 0.005);
    axBottom = axes(fig, "Units", "normalized", ...
        "Position", [0.14 0.10 0.80 0.36]);
    draw_panel_c(axBottom, panel, style, false, false, 0.01);
    exportgraphics(fig, path, "ContentType", "vector", ...
        "BackgroundColor", "white");
    close(fig);
end


function export_refined_supplementary(panelB, fullGroups, validGroups, ...
        panelC, style, outDir)
    export_si_local_full_ci(panelB, style, ...
        fullfile(outDir, "FigS_local_threshold_full_CI.pdf"));
    export_si_global_table(fullGroups, style, ...
        fullfile(outDir, "FigS_global_threshold_all_groups.pdf"), ...
        "All formal groups", true);
    export_si_global_table(panelC.weaklyIdentified, style, ...
        fullfile(outDir, "FigS_global_threshold_weakly_identified.pdf"), ...
        "Weakly identified PGG-valid groups", false);
    export_si_coefficients(validGroups, style, ...
        fullfile(outDir, "FigS_coefficient_forest_full.pdf"));
end


function export_si_local_full_ci(panel, style, path)
    fig = standalone_figure([18.0, 12.0]);
    layout = tiledlayout(fig, 2, 2, "Padding", "compact", ...
        "TileSpacing", "compact");
    for index = 1:numel(panel.displayR)
        ax = nexttile(layout, index); hold(ax, "on"); box(ax, "on");
        R = panel.displayR(index);
        T = panel.data(abs(panel.data.R - R) < 1e-10, :);
        theoryX = linspace(0.18, 0.82, 100);
        plot(ax, theoryX, theoryX * (1 - R), "-", ...
            "Color", style.rColors(index, :), "LineWidth", 1.6);
        errorbar(ax, T.x, T.estimate, T.estimate - T.lower_CI, ...
            T.upper_CI - T.estimate, "LineStyle", "none", ...
            "Marker", style.rMarkers(index), "MarkerSize", 4.5, ...
            "MarkerFaceColor", style.rColors(index, :), ...
            "MarkerEdgeColor", [0.15 0.15 0.15], ...
            "Color", style.rColors(index, :), "CapSize", 2);
        title(ax, sprintf("R=%.2f", R), "FontWeight", "normal");
        xlabel(ax, "Cooperator fraction, x"); ylabel(ax, "Local threshold");
        xlim(ax, [0.18 0.82]); apply_style(ax, style);
    end
    exportgraphics(fig, path, "ContentType", "vector", ...
        "BackgroundColor", "white"); close(fig);
end


function export_si_global_table(T, style, path, titleText, markInvalid)
    fig = standalone_figure([9.2, 7.2]);
    ax = axes(fig); hold(ax, "on"); box(ax, "on");
    if isempty(T)
        text(ax, 0.5, 0.5, "No groups met the predefined weak-identification criteria", ...
            "Units", "normalized", "HorizontalAlignment", "center");
        axis(ax, "off");
    else
        theory = 1 - T.R;
        limits = data_limits([theory; T.ci_low; T.ci_high], 0.06);
        plot(ax, limits, limits, "k--", "LineWidth", 1.8);
        games = sortrows(unique(T(:, ["k","n","M"]), "rows"), ...
            ["k","n","M"]);
        epsilonValues = sort(unique(T.epsilon));
        for row = 1:height(T)
            gameIndex = find(T.k(row) == games.k & T.n(row) == games.n & ...
                T.M(row) == games.M, 1);
            epsilonIndex = find(abs(T.epsilon(row) - epsilonValues) < 1e-12, 1);
            if isscalar(epsilonValues)
                epsilonIndex = 1 + (T.epsilon(row) > 0.0075);
            end
            errorbar(ax, theory(row), T.P_star_hat_x1(row), ...
                T.P_star_hat_x1(row) - T.ci_low(row), ...
                T.ci_high(row) - T.P_star_hat_x1(row), ...
                "LineStyle", "none", "Marker", style.gameMarkers(gameIndex), ...
                "MarkerSize", 5, "MarkerFaceColor", ...
                style.epsilonColors(epsilonIndex, :), ...
                "MarkerEdgeColor", [0.12 0.12 0.12], ...
                "Color", style.epsilonColors(epsilonIndex, :), "CapSize", 2);
        end
        if markInvalid && ismember("isPGGValid", string(T.Properties.VariableNames))
            invalid = ~T.isPGGValid;
            plot(ax, theory(invalid), T.P_star_hat_x1(invalid), "ko", ...
                "MarkerSize", 8, "LineWidth", 1.0, ...
                "MarkerFaceColor", "none", "HandleVisibility", "off");
        end
        xlabel(ax, "Theoretical threshold, 1-R");
        ylabel(ax, "Estimated global threshold");
        xlim(ax, limits); ylim(ax, limits); apply_style(ax, style);
    end
    title(ax, titleText, "FontWeight", "normal");
    exportgraphics(fig, path, "ContentType", "vector", ...
        "BackgroundColor", "white"); close(fig);
end


function paths = export_standalone_panels(panelA, panelB, panelC, panelD, style, outDir)
    paths = [
        string(fullfile(outDir, "Fig2a_equation_collapse.pdf"))
        string(fullfile(outDir, "Fig2b_local_boundary.pdf"))
        string(fullfile(outDir, "Fig2c_global_threshold.pdf"))
        string(fullfile(outDir, "Fig2d_coefficient_validation.pdf"))
    ];
    fig = standalone_figure([9.2, 7.2]);
    ax = axes(fig); draw_panel_a(ax, panelA, style, true);
    exportgraphics(fig, paths(1), "ContentType", "vector", "BackgroundColor", "white");
    close(fig); drawnow;

    fig = standalone_figure([9.2, 7.2]);
    ax = axes(fig); draw_panel_b(ax, panelB, style, true);
    exportgraphics(fig, paths(2), "ContentType", "vector", "BackgroundColor", "white");
    close(fig); drawnow;

    fig = standalone_figure([9.2, 7.2]);
    ax = axes(fig); draw_panel_c(ax, panelC, style, true, false);
    exportgraphics(fig, paths(3), "ContentType", "vector", "BackgroundColor", "white");
    close(fig); drawnow;

    fig = standalone_figure([9.2, 7.2]);
    ax = axes(fig); draw_panel_d(ax, panelD, style, true);
    exportgraphics(fig, paths(4), "ContentType", "vector", "BackgroundColor", "white");
    close(fig); drawnow;
end


function paths = export_supplementary_figures(validCells, fullCells, panelB, ...
        fullGroups, validGroups, balancedGroups, style, outDir, forceExport, selected)
    if nargin < 9, forceExport = false; end
    if nargin < 10 || isempty(selected), selected = 1:8; end
    selected = unique(selected(:).');
    paths = [
        string(fullfile(outDir, "FigS_sparse_equation_all_cells.pdf"))
        string(fullfile(outDir, "FigS_sparse_equation_by_R.pdf"))
        string(fullfile(outDir, "FigS_sparse_local_boundary_small_multiples.pdf"))
        string(fullfile(outDir, "FigS_sparse_RMSE_by_R.pdf"))
        string(fullfile(outDir, "FigS_sparse_full_vs_balanced.pdf"))
        string(fullfile(outDir, "FigS_sparse_all_groups.pdf"))
        string(fullfile(outDir, "FigS_sparse_coefficients_by_R.pdf"))
        string(fullfile(outDir, "FigS_sparse_phase_maps.pdf"))
    ];
    labels = ["all-cell equation collapse", "equation collapse by R", ...
        "local-boundary small multiples", "RMSE by R", ...
        "full versus balanced", "all groups including r<=1", ...
        "coefficients by R", "phase maps"];
    exporters = {@() export_si_equation_all(validCells, style, paths(1)), ...
        @() export_si_equation_by_r(validCells, style, paths(2)), ...
        @() export_si_local_small(panelB, style, paths(3)), ...
        @() export_si_rmse(validGroups, balancedGroups, style, paths(4)), ...
        @() export_si_full_balanced(validGroups, balancedGroups, style, paths(5)), ...
        @() export_si_all_groups(fullGroups, style, paths(6)), ...
        @() export_si_coefficients(validGroups, style, paths(7)), ...
        @() export_si_phase_maps(fullCells, style, paths(8))};
    for index = selected
        fprintf("Exporting legacy SI %d/8: %s\n", index, labels(index));
        if forceExport || ~is_nonempty_file(paths(index))
            exporters{index}();
        end
    end
end


function export_si_equation_all(cells, style, path)
    fig = standalone_figure([9.2, 7.2]);
    ax = axes(fig); hold(ax, "on"); box(ax, "on");
    Z = cells.mean_P_actual - cells.mean_x_after_relax .* (1 - cells.R);
    Y = cells.mean_G_over_epsilon;
    limits = data_limits([Z; Y], 0.05);
    plot(ax, limits, limits, "k--", "LineWidth", style.theoryWidth);
    scatter(ax, Z, Y, 10, cells.R, "filled", "MarkerFaceAlpha", 0.35);
    colorbar(ax); colormap(ax, parula(128));
    xlabel(ax, "$Z=P-x(1-R)$", "Interpreter", "latex");
    ylabel(ax, "$\widehat G/\epsilon$", "Interpreter", "latex");
    title(ax, "All PGG-valid cell means", "FontWeight", "normal");
    xlim(ax, limits); ylim(ax, limits); axis(ax, "square"); apply_style(ax, style);
    exportgraphics(fig, path, "ContentType", "vector");
    save_legacy_supplementary_fig(fig, path, style); close(fig);
end


function export_si_equation_by_r(cells, style, path)
    rValues = sort(unique(cells.R));
    colorPosition = round(linspace(1, size(style.rColormap, 1), numel(rValues)));
    colors = style.rColormap(colorPosition, :);
    fig = standalone_figure([18.5, 12.0]);
    layout = tiledlayout(fig, 2, 3, "Padding", "compact", "TileSpacing", "compact");
    for index = 1:numel(rValues)
        ax = nexttile(layout, index); hold(ax, "on"); box(ax, "on");
        T = cells(abs(cells.R - rValues(index)) < 1e-10, :);
        Z = T.mean_P_actual - T.mean_x_after_relax .* (1 - T.R);
        Y = T.mean_G_over_epsilon; limits = data_limits([Z; Y], 0.05);
        plot(ax, limits, limits, "k--", "LineWidth", style.theoryWidth);
        scatter(ax, Z, Y, 11, colors(index, :), "filled", ...
            "MarkerFaceAlpha", 0.38, "MarkerEdgeColor", "none");
        title(ax, sprintf("R=%.2f", rValues(index)), "FontWeight", "normal");
        xlabel(ax, "Z"); ylabel(ax, "G/epsilon");
        xlim(ax, limits); ylim(ax, limits); axis(ax, "square"); apply_style(ax, style);
    end
    exportgraphics(fig, path, "ContentType", "vector");
    save_legacy_supplementary_fig(fig, path, style); close(fig);
end


function export_si_local_small(panel, style, path)
    % This legacy SI contains every fitted R, whereas the current main panel
    % deliberately displays only four R values.  Derive the layout, colors,
    % and markers from the actual SI data rather than indexing the four-item
    % main-panel style arrays.
    rValues = panel.rValues(:).';
    numberOfR = numel(rValues);
    numberOfColumns = min(3, numberOfR);
    numberOfRows = ceil(numberOfR / numberOfColumns);
    fig = standalone_figure([18.5, 5.2 * numberOfRows]);
    layout = tiledlayout(fig, numberOfRows, numberOfColumns, ...
        "Padding", "compact", "TileSpacing", "compact");
    colorPosition = round(linspace(1, size(style.rColormap, 1), numberOfR));
    colors = style.rColormap(colorPosition, :);
    markers = ["o", "s", "^", "d", "v", ">", "<"];
    for index = 1:numberOfR
        ax = nexttile(layout, index); hold(ax, "on"); box(ax, "on");
        R = rValues(index); T = panel.data(abs(panel.data.R - R) < 1e-10, :);
        if isempty(T)
            warning("SparseFinalFig2:EmptyLegacyLocalR", ...
                "No local-threshold rows found for R=%.3f.", R);
            axis(ax, "off");
            continue;
        end
        xTheory = linspace(min(T.x), max(T.x), 100);
        plot(ax, xTheory, xTheory * (1 - R), "-", ...
            "Color", colors(index, :), "LineWidth", style.theoryWidth);
        errorbar(ax, T.x, T.estimate, T.estimate - T.lower_CI, ...
            T.upper_CI - T.estimate, "LineStyle", "none", ...
            "Marker", markers(mod(index - 1, numel(markers)) + 1), ...
            "MarkerSize", 4.5, "MarkerFaceColor", colors(index, :), ...
            "MarkerEdgeColor", [0.15 0.15 0.15], ...
            "Color", colors(index, :), "CapSize", 2);
        title(ax, sprintf("R=%.2f", R), "FontWeight", "normal");
        xlabel(ax, "x"); if index == 1, ylabel(ax, "Local threshold"); end
        apply_style(ax, style);
    end
    exportgraphics(fig, path, "ContentType", "vector");
    save_legacy_supplementary_fig(fig, path, style); close(fig);
end


function export_si_rmse(full, balanced, style, path)
    fig = standalone_figure([9.2, 7.2]); ax = axes(fig); hold(ax, "on"); box(ax, "on");
    rValues = sort(unique(full.R)); fullRMSE = nan(size(rValues)); balRMSE = fullRMSE;
    for index = 1:numel(rValues)
        F = full(abs(full.R - rValues(index)) < 1e-10, :);
        B = balanced(abs(balanced.R - rValues(index)) < 1e-10 & balanced.isPGGValid, :);
        fullRMSE(index) = sqrt(mean((F.P_star_hat_x1 - (1 - F.R)).^2));
        balRMSE(index) = sqrt(mean((B.P_star_hat_x1 - (1 - B.R)).^2));
    end
    plot(ax, rValues, fullRMSE, "-o", "LineWidth", 1.5, ...
        "MarkerFaceColor", style.epsilonColors(1, :));
    plot(ax, rValues, balRMSE, "--s", "LineWidth", 1.4, ...
        "MarkerFaceColor", style.epsilonColors(2, :));
    xlabel(ax, "R"); ylabel(ax, "Global-threshold RMSE");
    legend(ax, {"full", "balanced300"}, "Box", "off"); apply_style(ax, style);
    exportgraphics(fig, path, "ContentType", "vector");
    save_legacy_supplementary_fig(fig, path, style); close(fig);
end


function export_si_full_balanced(full, balanced, style, path)
    B = balanced(balanced.isPGGValid, :);
    F = sortrows(full, ["R", "epsilon", "k", "n", "M"]);
    B = sortrows(B, ["R", "epsilon", "k", "n", "M"]);
    fig = standalone_figure([9.2, 7.2]); ax = axes(fig); hold(ax, "on"); box(ax, "on");
    limits = data_limits([F.P_star_hat_x1; B.P_star_hat_x1], 0.06);
    plot(ax, limits, limits, "k--", "LineWidth", style.theoryWidth);
    rValues = sort(unique(F.R));
    colors = lines(numel(rValues));
    handles = gobjects(numel(rValues), 1);
    for index = 1:numel(rValues)
        mask = abs(F.R - rValues(index)) < 1e-10;
        handles(index) = plot(ax, F.P_star_hat_x1(mask), ...
            B.P_star_hat_x1(mask), "o", "LineStyle", "none", ...
            "MarkerSize", 4.5, "MarkerFaceColor", colors(index, :), ...
            "MarkerEdgeColor", [0.15 0.15 0.15]);
    end
    legend(ax, handles, cellstr(compose("R=%.2g", rValues)), ...
        "Location", "northwest", "NumColumns", 2, "Box", "off", ...
        "FontSize", style.smallFont);
    xlabel(ax, "Full-data estimate"); ylabel(ax, "Balanced300 estimate");
    xlim(ax, limits); ylim(ax, limits); axis(ax, "square"); apply_style(ax, style);
    exportgraphics(fig, path, "ContentType", "vector");
    save_legacy_supplementary_fig(fig, path, style); close(fig);
end


function export_si_all_groups(full, style, path)
    T = full; T.theory = 1 - T.R; T.covered = T.ci_low <= T.theory & T.ci_high >= T.theory;
    % draw_panel_c now reads displayData because the main-text panel can
    % distinguish the complete table from its well-identified subset.  This
    % legacy SI intentionally displays the complete table, including r<=1,
    % so both fields point to T here.
    panel.data = T;
    panel.displayData = T;
    panel.groups = height(T); panel.coverage = sum(T.covered);
    panel.rmse = sqrt(mean((T.P_star_hat_x1 - T.theory).^2));
    fig = standalone_figure([9.2, 7.2]); ax = axes(fig);
    draw_panel_c(ax, panel, style, false, true);
    invalid = ~T.isPGGValid;
    hold(ax, "on"); plot(ax, T.theory(invalid), T.P_star_hat_x1(invalid), ...
        "rx", "MarkerSize", 8, "LineWidth", 1.4);
    text(ax, 0.97, 0.05, sprintf("all groups %d; invalid %d", ...
        height(T), sum(invalid)), "Units", "normalized", ...
        "HorizontalAlignment", "right", "FontSize", style.smallFont);
    exportgraphics(fig, path, "ContentType", "vector");
    save_legacy_supplementary_fig(fig, path, style); close(fig);
end


function export_si_coefficients(T, style, path)
    T.delta_alpha = T.intercept_over_epsilon;
    T.delta_b = T.coef_p - 1;
    T.delta_s = T.state_slope_hat - (1 - T.R);
    lows = {T.intercept_ci_low, T.coef_p_ci_low - 1, ...
        T.state_slope_ci_low - (1 - T.R)};
    highs = {T.intercept_ci_high, T.coef_p_ci_high - 1, ...
        T.state_slope_ci_high - (1 - T.R)};
    estimates = {T.delta_alpha, T.delta_b, T.delta_s};
    labels = {"alpha", "b-1", "s-(1-R)"};
    fig = standalone_figure([18.0, 6.0]);
    layout = tiledlayout(fig, 1, 3, "Padding", "compact", "TileSpacing", "compact");
    for category = 1:3
        ax = nexttile(layout, category); hold(ax, "on"); box(ax, "on"); yline(ax, 0, "k--");
        errorbar(ax, T.R, estimates{category}, estimates{category} - lows{category}, ...
            highs{category} - estimates{category}, "LineStyle", "none", ...
            "Marker", "o", "MarkerSize", 3, "Color", [0.35 0.35 0.35], ...
            "CapSize", 1.5);
        xlabel(ax, "R"); ylabel(ax, labels{category}); apply_style(ax, style);
    end
    exportgraphics(fig, path, "ContentType", "vector");
    save_legacy_supplementary_fig(fig, path, style); close(fig);
end


function export_si_phase_maps(cells, style, path)
    rValues = [0.60, 0.70, 0.80, 0.85, 0.90];
    T = cells(cells.analysis == "full" & cells.k == 6 & cells.n == 3 & ...
        cells.M == 15 & abs(cells.epsilon - 0.01) < 1e-12 & ...
        ismembertol(cells.R, rValues, 1e-10), :);
    maxAbs = max(abs(T.mean_G_over_epsilon), [], "omitnan");
    if ~isfinite(maxAbs) || maxAbs <= 0, maxAbs = 1; end
    negativeColor = hex2rgb('#3A5F8A');
    positiveColor = hex2rgb('#E4985D');
    continuousMap = diverging_colormap(negativeColor, [1 1 1], ...
        positiveColor, 257);
    binaryMap = [negativeColor; positiveColor];

    % The sixth column is reserved for one shared colorbar per row.  Each
    % data axis retains its own colormap, so the continuous and binary maps
    % cannot be confused when the FIG is edited interactively.
    fig = standalone_figure([19.0, 8.8]);
    layout = tiledlayout(fig, 2, 6, "Padding", "compact", ...
        "TileSpacing", "compact");
    topAxes = gobjects(5, 1);
    bottomAxes = gobjects(5, 1);
    for index = 1:5
        group = T(abs(T.R - rValues(index)) < 1e-10, :);
        topAxes(index) = nexttile(layout, index);
        ax = topAxes(index); hold(ax, "on"); box(ax, "on");
        scatter(ax, group.mean_x_after_relax, group.mean_P_actual, 42, ...
            group.mean_G_over_epsilon, "s", "filled", ...
            "MarkerEdgeColor", [0.82 0.82 0.82], "LineWidth", 0.25, ...
            "HandleVisibility", "off");
        colormap(ax, continuousMap); clim(ax, [-maxAbs maxAbs]);
        theoryX = linspace(min(group.mean_x_after_relax), ...
            max(group.mean_x_after_relax), 120);
        plot(ax, theoryX, theoryX .* (1 - rValues(index)), "k-", ...
            "LineWidth", 1.25, "DisplayName", "$P^*(x)=x(1-R)$");
        title(ax, sprintf("R=%.2f", rValues(index)), "FontWeight", "normal");
        if index == 1
            ylabel(ax, "Normalized steering, $P$", "Interpreter", "latex");
            legend(ax, "Location", "northwest", "Box", "off", ...
                "Interpreter", "latex", "FontSize", style.smallFont);
        end
        ax.XTickLabel = [];
        apply_style(ax, style);

        bottomAxes(index) = nexttile(layout, 6 + index);
        ax = bottomAxes(index); hold(ax, "on"); box(ax, "on");
        scatter(ax, group.mean_x_after_relax, group.mean_P_actual, 42, ...
            double(group.mean_G_over_epsilon >= 0), "s", "filled", ...
            "MarkerEdgeColor", [0.82 0.82 0.82], "LineWidth", 0.25, ...
            "HandleVisibility", "off");
        colormap(ax, binaryMap); clim(ax, [-0.5 1.5]);
        plot(ax, theoryX, theoryX .* (1 - rValues(index)), "k-", ...
            "LineWidth", 1.25, "HandleVisibility", "off");
        xlabel(ax, "Relaxed cooperator fraction, $x$", ...
            "Interpreter", "latex");
        if index == 1
            ylabel(ax, "Normalized steering, $P$", "Interpreter", "latex");
        end
        apply_style(ax, style);
    end

    continuousBar = colorbar(topAxes(end));
    continuousBar.Layout.Tile = 6;
    continuousBar.Label.String = "$\widehat{G}/\epsilon$";
    continuousBar.Label.Interpreter = "latex";
    continuousBar.FontName = style.font;
    continuousBar.FontSize = style.smallFont;

    binaryBar = colorbar(bottomAxes(end));
    binaryBar.Layout.Tile = 12;
    binaryBar.Ticks = [0 1];
    binaryBar.TickLabels = {"Negative drift", "Non-negative drift"};
    binaryBar.FontName = style.font;
    binaryBar.FontSize = style.smallFont;

    exportgraphics(fig, path, "ContentType", "vector");
    save_legacy_supplementary_fig(fig, path, style); close(fig);
end


function cmap = diverging_colormap(lowColor, centerColor, highColor, numberOfColors)
    numberOfColors = max(3, round(numberOfColors));
    anchors = [0 0.5 1];
    cmap = interp1(anchors, [lowColor; centerColor; highColor], ...
        linspace(0, 1, numberOfColors), "linear");
    cmap = min(max(cmap, 0), 1);
end


function save_legacy_supplementary_fig(fig, pdfPath, style)
    if ~style.cfg.export_legacy_supplementary_fig
        return;
    end
    [folder, name] = fileparts(char(pdfPath));
    figPath = fullfile(folder, [name '.fig']);
    % standalone_figure is intentionally invisible during batch export.
    % FIG files preserve the Visible property, so saving that state makes a
    % valid file appear not to open when double-clicked.  Store the editable
    % copy as visible, then restore the batch figure state before closing it.
    previousVisibility = fig.Visible;
    fig.Visible = 'on';
    drawnow;
    savefig(fig, figPath);
    fig.Visible = previousVisibility;
    fprintf("  editable FIG: %s\n", figPath);
end


function fig = standalone_figure(sizeCm)
    fig = figure("Color", "white", "Units", "centimeters", ...
        "Position", [1, 1, sizeCm(1), sizeCm(2)], ...
        "PaperPositionMode", "auto", "Visible", "off");
end


function ready = is_nonempty_file(path)
    if ~isfile(path)
        ready = false;
        return;
    end
    information = dir(path);
    ready = ~isempty(information) && information(1).bytes > 0;
end


function limits = data_limits(values, fraction)
    values = values(isfinite(values));
    width = max(values) - min(values);
    if width <= 0, width = max(1, abs(values(1))); end
    padding = fraction * width;
    limits = [min(values) - padding, max(values) + padding];
end


function apply_style(ax, style)
    ax.FontName = style.font; ax.FontSize = style.axisFont;
    ax.LineWidth = style.axisWidth; ax.TickDir = "out";
    ax.TickLength = [0.018 0.018]; ax.Layer = "top"; box(ax, "on");
    ax.XLabel.FontSize = style.labelFont; ax.YLabel.FontSize = style.labelFont;
end


function add_panel_label(ax, label, style)
    text(ax, -0.12, 1.12, label, "Units", "normalized", ...
        "HorizontalAlignment", "left", "VerticalAlignment", "top", ...
        "FontName", style.font, "FontSize", style.panelFont, ...
        "FontWeight", "bold", "Clipping", "off");
end


function T = empty_plot_table()
    T = table('Size', [0, 15], ...
        'VariableTypes', ["string", repmat("double", 1, 11), "string", "string", "string"], ...
        'VariableNames', ["panel", "k", "n", "M", "R", "epsilon", ...
        "x", "P", "estimate", "lower_CI", "upper_CI", "theory", ...
        "point_type", "source_file", "notes"]);
end


function T = append_plot_row(T, panel, k, n, M, R, epsilon, x, P, ...
        estimate, lowerCI, upperCI, theory, pointType, sourceFile)
    row = table(string(panel), k, n, M, R, epsilon, x, P, estimate, ...
        lowerCI, upperCI, theory, string(pointType), string(sourceFile), "", ...
        'VariableNames', T.Properties.VariableNames);
    T = [T; row];
end


function write_statistics(path, panelA, panelB, panelC, panelD, ...
        fullGroups, validGroups, excludedGroups, paths, gAudit)
    fileId = fopen(path, "w");
    if fileId < 0, error("SparseFinalFig2:StatisticsWrite", "Cannot write %s", path); end
    closer = onCleanup(@() fclose(fileId)); %#ok<NASGU>
    fprintf(fileId, "Final sparse-steering Figure 2 statistics\n\n");
    fprintf(fileId, "Selection rule: r0=k*n/(k+n-1); r=R*r0; main text requires r>1+1e-12.\n");
    fprintf(fileId, "No R2, RMSE, coverage, or theory-agreement filtering.\n");
    fprintf(fileId, "Full pooled estimator is primary; balanced300 appears only in Supplementary.\n\n");
    fprintf(fileId, "G_hat normalization audit: passed=%d\nformula=%s\nrunner=%s\n\n", ...
        gAudit.correctStructureFactor, strjoin(gAudit.formula, ""), gAudit.runnerPath);
    inventory = readtable(paths.inventory, "TextType", "string");
    inventory = inventory(inventory.analysis == "full" & inventory.r > 1 + 1e-12, :);
    fprintf(fileId, "Computation inventory (PGG-valid full design)\n");
    fprintf(fileId, "total_simulation_runs=%d\naccepted_trajectories=%d\ncell_estimates=%d\n\n", ...
        sum(inventory.num_launched), sum(inventory.num_accepted), panelA.cells);
    fprintf(fileId, "Panel a\n");
    fprintf(fileId, "total_cells=%d\naccepted_trajectories=%d\n", panelA.cells, panelA.accepted);
    fprintf(fileId, "alpha=%.8f CI=[%.8f,%.8f]\n", panelA.alpha, panelA.alphaCI);
    fprintf(fileId, "beta=%.8f CI=[%.8f,%.8f]\n", panelA.beta, panelA.betaCI);
    fprintf(fileId, "weighted_R2=%.8f\nweighted_RMSE=%.8f\n", ...
        panelA.weightedR2, panelA.weightedRMSE);
    fprintf(fileId, "binning=%d equal-cell-count Z bins; inverse-variance weighted means; cell bootstrap CI\n", ...
        panelA.numberOfBins);
    fprintf(fileId, "number_of_bins=%d\n\n", panelA.numberOfBins);
    fprintf(fileId, "Panel b\nR values: %s\n", mat2str(panelB.rValues));
    fprintf(fileId, "main_figure_R_values: %s\n", mat2str(panelB.displayR));
    fprintf(fileId, "Plotting note: intervals extending beyond the plotting range are indicated by arrows; numerical CIs are unchanged.\n");
    fprintf(fileId, "Caption text: Intervals extending beyond the plotting range are indicated by arrows; full bootstrap intervals are reported in Supplementary Fig. Sx.\n");
    fprintf(fileId, "Weak-identification criteria: slope CI contains zero; crossing outside sampled actual-P range; finite bootstrap fraction below 0.80; nonfinite CI; or fewer than 3 actual-P cells.\n");
    fprintf(fileId, "x_coordinate=accepted-count weighted mean_x_after_relax\n");
    fprintf(fileId, "bootstrap=%s\n", panelB.bootstrapMethod);
    fprintf(fileId, "weakly_identified_main_R_points=%d\nclipped_main_CIs=%d\n", ...
        height(panelB.weaklyIdentified), panelB.clippedCount);
    for row = 1:height(panelB.weaklyIdentified)
        W = panelB.weaklyIdentified;
        fprintf(fileId, "  weak local: R=%.2f x0=%.2f estimate=%.8f reason=%s\n", ...
            W.R(row), W.x0_target(row), W.estimate(row), W.weak_reason(row));
    end
    for R = panelB.rValues
        group = panelB.data(abs(panelB.data.R - R) < 1e-10, :);
        errors = group.estimate - group.theory;
        fprintf(fileId, "R=%.2f thresholds=%d median_R2=%.8f RMSE=%.8f MAE=%.8f coverage=%d/%d\n", ...
            R, height(group), median(group.local_fit_R2), sqrt(mean(errors.^2)), ...
            mean(abs(errors)), sum(group.covered), height(group));
    end
    fprintf(fileId, "\nPanel c\nformal_groups_total=%d\nPGG_valid=%d\nexcluded=%d\n", ...
        height(fullGroups), height(validGroups), height(excludedGroups));
    fprintf(fileId, "coverage=%d/%d\nRMSE=%.8f\nMAE=%.8f\n", ...
        panelC.coverage, panelC.groups, panelC.rmse, panelC.mae);
    fprintf(fileId, "well_identified_for_main_figure=%d; weakly_identified_for_supplementary=%d\n", ...
        height(panelC.displayData), height(panelC.weaklyIdentified));
    fprintf(fileId, "Well-identification criteria are independent of theory agreement and use the same slope/sampled-range/bootstrap-finiteness logic as Panel b.\n");
    fprintf(fileId, "Plotting note: intervals extending beyond the plotting range are indicated by arrows; numerical CIs are unchanged.\n");
    fprintf(fileId, "uncovered_groups:\n");
    for row = 1:height(panelC.uncovered)
        T = panelC.uncovered;
        fprintf(fileId, "  k=%d n=%d M=%d R=%.2f epsilon=%.3g estimate=%.8f CI=[%.8f,%.8f] theory=%.8f\n", ...
            T.k(row), T.n(row), T.M(row), T.R(row), T.epsilon(row), ...
            T.P_star_hat_x1(row), T.ci_low(row), T.ci_high(row), T.theory(row));
    end
    fprintf(fileId, "\nExcluded groups:\n");
    fprintf(fileId, "%s\n", evalc("disp(excludedGroups)"));
    fprintf(fileId, "\nPanel d\nalpha_coverage=%d/%d\nb_coverage=%d/%d\ns_coverage=%d/%d\n", ...
        panelD.alphaCoverage, panelD.groups, panelD.bCoverage, panelD.groups, ...
        panelD.sCoverage, panelD.groups);
    fprintf(fileId, "pooled_deviation_summaries (IVW mean and normal 95%% CI):\n");
    fprintf(fileId, "%s\n", evalc("disp(panelD.summary)"));
    fprintf(fileId, "\nInput files:\n%s\n%s\n%s\n", ...
        paths.cells, paths.bootstrap, paths.inventory);
    if strlength(paths.checklist) > 0, fprintf(fileId, "%s\n", paths.checklist); end
end
