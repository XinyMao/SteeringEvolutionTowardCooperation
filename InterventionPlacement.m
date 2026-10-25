function fig = plot_figure3d_exposure_p003(cfg)
%PLOT_FIGURE3D_EXPOSURE_P003 Editable Figure 3d using the final p=0.03 data.
%
%   fig = plot_figure3d_exposure_p003()
%   fig = plot_figure3d_exposure_p003(cfg)
%
% The main inferential unit is a configuration mean (30 trajectories). Raw
% simulations are never modified. The solid fit uses all configurations; the
% grey dashed fit excludes clustered placement as a pre-specified sensitivity
% analysis. Protocol means and their two-dimensional 95% intervals are read
% from the completed network-bootstrap summary.

if nargin < 1 || isempty(cfg), cfg = struct(); end
cfg = defaults(cfg);

thisFile = mfilename('fullpath');
defaultRoot = fileparts(fileparts(fileparts(thisFile)));
if strlength(string(cfg.project_root)) == 0
    rootDir = defaultRoot;
else
    rootDir = char(cfg.project_root);
end
resultDir = fullfile(rootDir,'results', ...
    'exp04_finite_exposure_mechanism_N5000_p003_w0003_tau5_case64');
figureDir = fullfile(rootDir,'figures');
if ~exist(figureDir,'dir'), mkdir(figureDir); end

configFile = fullfile(resultDir,'finite_exposure_p003_configurations.csv');
summaryFile = fullfile(resultDir,'finite_exposure_p003_protocol_summary.csv');
statsFile = fullfile(resultDir,'finite_exposure_p003_statistics.json');
required = {configFile,summaryFile,statsFile};
for i=1:numel(required)
    if ~isfile(required{i}), error('Missing final p=0.03 data: %s',required{i}); end
end

D = readtable(configFile,'TextType','string');
S = readtable(summaryFile,'TextType','string');
J = jsondecode(fileread(statsFile));
validateData(D,S,J);

pal.blue = hex2rgb('#3A5F8A');
pal.lightBlue = hex2rgb('#79A7C8');
pal.green = hex2rgb('#63A69F');
pal.orange = hex2rgb('#E4985D');
pal.magenta = hex2rgb('#A65C85');
pal.grey = hex2rgb('#B8BDC6');
pal.dark = [0.08 0.08 0.08];

protocols = ["random","d_interface","metropolis_exposure","c_core","clustered"];
labels = {'Random','D-interface','Exposure-biased','C-core','Clustered'};
colors = [pal.grey;pal.blue;pal.green;pal.orange;pal.magenta];
markers = {'o','^','d','s','v'};

% The JSON values are the final audited all-configuration statistics.
allFit.intercept = J.intercept;
allFit.slope = J.slope;
allFit.low = J.slope_CI_low;
allFit.high = J.slope_CI_high;
allFit.R2 = J.R2;

nonCluster = D(D.placement_protocol~="clustered",:);
nonFit = clusteredRegression(nonCluster,cfg.bootstrap,cfg.seed);

set(groot,'defaultTextInterpreter','latex', ...
    'defaultAxesTickLabelInterpreter','latex', ...
    'defaultLegendInterpreter','latex');
fig = figure('Color','w','Units','centimeters','Position',[2 2 10.2 8.2], ...
    'Name','Figure 3d: exposure mechanism (p=0.03)','NumberTitle','off');
% Independent axes (not tiledlayout) remain directly resizable in Property
% Inspector. Positions are normalized and can be edited after opening the FIG.
ax = axes('Parent',fig,'Units','normalized','Position',[0.14 0.14 0.72 0.78]);
hold(ax,'on');

states = sort(unique(D.x_target));
stateMap = interp1([min(states),median(states),max(states)], ...
    [pal.blue;pal.orange;pal.green],linspace(min(states),max(states),256),'linear');
pointColors = interp1(linspace(min(states),max(states),size(stateMap,1)), ...
    stateMap,D.x_target,'linear');
hCells = scatter(ax,D.delta_eta0,D.delta_drift,cfg.point_size,pointColors,'o','filled', ...
    'MarkerFaceAlpha',cfg.point_alpha,'MarkerEdgeColor','none', ...
    'DisplayName','Configuration estimates');

xAll = linspace(min(D.delta_eta0),max(D.delta_eta0),240);
hAll = plot(ax,xAll,allFit.intercept+allFit.slope*xAll,'k-', ...
    'LineWidth',1.65,'DisplayName','All-configuration fit');
xNC = linspace(min(nonCluster.delta_eta0),max(nonCluster.delta_eta0),180);
hNC = plot(ax,xNC,nonFit.intercept+nonFit.slope*xNC,'--', ...
    'Color',[.42 .42 .42],'LineWidth',1.25, ...
    'DisplayName','Fit excluding Clustered');

meanHandles = gobjects(numel(protocols),1);
for ip=1:numel(protocols)
    q=S(S.placement_protocol==protocols(ip),:);
    if height(q)~=1, error('Expected one summary row for %s.',protocols(ip)); end
    drawXYInterval(ax,q.mean_delta_eta0,q.mean_delta_drift, ...
        q.eta_CI_low,q.eta_CI_high,q.drift_CI_low,q.drift_CI_high, ...
        colors(ip,:),markers{ip},7.2,pal.dark);
    meanHandles(ip)=plot(ax,nan,nan,markers{ip},'LineStyle','none', ...
        'MarkerSize',6.6,'MarkerFaceColor',colors(ip,:), ...
        'MarkerEdgeColor',pal.dark,'LineWidth',.8,'DisplayName',labels{ip});
end

xline(ax,0,':','Color',[.65 .65 .65],'LineWidth',.85,'HandleVisibility','off');
yline(ax,0,':','Color',[.65 .65 .65],'LineWidth',.85,'HandleVisibility','off');
xlim(ax,[-0.00635 0.00165]); ylim(ax,[-0.0060 0.00245]);
xlabel(ax,'Exposure difference, $\Delta\eta_0$','FontSize',9.5);
ylabel(ax,'Drift difference, $\Delta\langle\dot{x}\rangle$','FontSize',9.5);
ax.FontName=cfg.font_name; ax.FontSize=8; ax.LineWidth=.9; ax.TickDir='out';
ax.Box='on'; ax.Layer='top'; grid(ax,'off');
ax.XAxis.Exponent=-3; ax.YAxis.Exponent=-3;

colormap(ax,stateMap); clim(ax,[min(states),max(states)]);
cb=colorbar(ax); cb.Position=[0.89 0.18 0.025 0.69];
cb.Ticks=states; cb.FontName=cfg.font_name; cb.FontSize=7;
cb.Label.String='Initial cooperation, $x_0$'; cb.Label.Interpreter='latex';

text(ax,-.12,1.04,'\textbf{d}','Units','normalized','FontSize',11, ...
    'Clipping','off');
betaAllText=['All: $\beta=' num2str(allFit.slope,'%.3f') ...
    '\,[' num2str(allFit.low,'%.3f') ',' num2str(allFit.high,'%.3f') ...
    ']$, $R^2=' num2str(allFit.R2,'%.3f') '$'];
betaNCText=['No-clust.: $\beta=' num2str(nonFit.slope,'%.3f') ...
    '\,[' num2str(nonFit.low,'%.3f') ',' num2str(nonFit.high,'%.3f') ...
    ']$, $R^2=' num2str(nonFit.R2,'%.3f') '$'];
text(ax,.98,.98,{'$p=0.03,\;r=2$',betaAllText,betaNCText}, ...
    'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top', ...
    'FontSize',5.2,'Color',[.18 .18 .18]);

legend(ax,[hCells;hAll;hNC;meanHandles],'Location','southwest', ...
    'NumColumns',2,'FontSize',5.7,'Box','off');

% Magnified protocol means excluding Clustered. This is a two-dimensional
% zoom, not a performance ranking; Clustered remains visible in the main axis.
inset = axes('Parent',fig,'Units','normalized','Position',[0.19 0.61 0.28 0.24]);
hold(inset,'on');
for ip=1:4
    q=S(S.placement_protocol==protocols(ip),:);
    drawXYInterval(inset,1e3*q.mean_delta_eta0,1e3*q.mean_delta_drift, ...
        1e3*q.eta_CI_low,1e3*q.eta_CI_high, ...
        1e3*q.drift_CI_low,1e3*q.drift_CI_high, ...
        colors(ip,:),markers{ip},5.5,pal.dark);
end
xline(inset,0,':','Color',[.67 .67 .67],'LineWidth',.65);
yline(inset,0,':','Color',[.67 .67 .67],'LineWidth',.65);
xlim(inset,[-.08 .08]); ylim(inset,[-.23 .12]);
xticks(inset,[-.05 0 .05]); yticks(inset,[-.2 -.1 0 .1]);
xlabel(inset,'$10^3\Delta\eta_0$','FontSize',5.5);
ylabel(inset,'$10^3\Delta\langle\dot{x}\rangle$','FontSize',5.5);
title(inset,'Non-clustered protocol means','FontWeight','normal','FontSize',5.8);
inset.FontName=cfg.font_name; inset.FontSize=5.4; inset.LineWidth=.65;
inset.TickDir='out'; inset.Box='on'; inset.Layer='top';

figPath=fullfile(figureDir,'Figure3d_exposure_p003_editable.fig');
pdfPath=fullfile(figureDir,'Figure3d_exposure_p003.pdf');
pngPath=fullfile(figureDir,'Figure3d_exposure_p003.png');
if cfg.export_fig, savefig(fig,figPath); end
if cfg.export_pdf, exportgraphics(fig,pdfPath,'ContentType','vector'); end
if cfg.export_png, exportgraphics(fig,pngPath,'Resolution',cfg.png_dpi); end

statsOut=fullfile(figureDir,'derived','fig3d_exposure_p003_plot_statistics.txt');
if ~exist(fileparts(statsOut),'dir'),mkdir(fileparts(statsOut));end
fid=fopen(statsOut,'w');
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'Figure 3d p=0.03 exposure mechanism\n');
fprintf(fid,'Configurations: %d; trajectories per configuration: %d\n',height(D),D.n_trajectories(1));
fprintf(fid,'All: slope %.9g, CI [%.9g, %.9g], R2 %.9g\n', ...
    allFit.slope,allFit.low,allFit.high,allFit.R2);
fprintf(fid,'Excluding Clustered: slope %.9g, CI [%.9g, %.9g], R2 %.9g\n', ...
    nonFit.slope,nonFit.low,nonFit.high,nonFit.R2);
fprintf(fid,'Source configurations: %s\nSource summaries: %s\n',configFile,summaryFile);

fprintf('\nEditable Figure 3d generated:\n  %s\n  %s\n  %s\n',figPath,pdfPath,pngPath);
fprintf('All configurations: slope %.4f, CI [%.4f, %.4f], R2 %.4f\n', ...
    allFit.slope,allFit.low,allFit.high,allFit.R2);
fprintf('Excluding Clustered: slope %.4f, CI [%.4f, %.4f], R2 %.4f\n', ...
    nonFit.slope,nonFit.low,nonFit.high,nonFit.R2);
end

function cfg=defaults(cfg)
cfg=setDefault(cfg,'project_root','');
cfg=setDefault(cfg,'font_name','Arial');
cfg=setDefault(cfg,'export_fig',true);
cfg=setDefault(cfg,'export_pdf',true);
cfg=setDefault(cfg,'export_png',true);
cfg=setDefault(cfg,'png_dpi',600);
cfg=setDefault(cfg,'point_alpha',0.075);
cfg=setDefault(cfg,'point_size',8);
cfg=setDefault(cfg,'bootstrap',1000);
cfg=setDefault(cfg,'seed',20260810);
end

function s=setDefault(s,name,value)
if ~isfield(s,name)||isempty(s.(name)),s.(name)=value;end
end

function validateData(D,S,J)
requiredD={'placement_protocol','x_target','network_seed','delta_eta0', ...
    'delta_drift','n_trajectories'};
requiredS={'placement_protocol','mean_delta_eta0','eta_CI_low','eta_CI_high', ...
    'mean_delta_drift','drift_CI_low','drift_CI_high'};
if ~all(ismember(requiredD,D.Properties.VariableNames))
    error('Configuration table columns: %s',strjoin(D.Properties.VariableNames,', '));
end
if ~all(ismember(requiredS,S.Properties.VariableNames))
    error('Summary table columns: %s',strjoin(S.Properties.VariableNames,', '));
end
if height(D)~=1408 || any(D.n_trajectories~=30)
    error('Incomplete p=0.03 configurations: rows=%d, trajectory counts=%s', ...
        height(D),mat2str(unique(D.n_trajectories)'));
end
if J.accepted_trajectories~=42240
    error('Expected 42,240 accepted trajectories; found %d.',J.accepted_trajectories);
end
end

function fit=clusteredRegression(T,B,seed)
x=T.delta_eta0; y=T.delta_drift;
X=[ones(height(T),1),x]; beta=X\y; fitted=X*beta;
fit.intercept=beta(1); fit.slope=beta(2);
fit.R2=1-sum((y-fitted).^2)/sum((y-mean(y)).^2);
networks=unique(T.network_seed); boot=nan(B,1);
rng(seed,'twister');
groups=cell(numel(networks),1);
for i=1:numel(networks),groups{i}=find(T.network_seed==networks(i));end
for b=1:B
    sampled=randi(numel(groups),numel(groups),1);
    rows=vertcat(groups{sampled});
    xb=[ones(numel(rows),1),x(rows)]; bb=xb\y(rows); boot(b)=bb(2);
end
ci=prctile(boot,[2.5 97.5]); fit.low=ci(1);fit.high=ci(2);
end

function drawXYInterval(ax,x,y,xLow,xHigh,yLow,yHigh,color,marker,size,edge)
plot(ax,[xLow xHigh],[y y],'-','Color',edge,'LineWidth',1.05,'HandleVisibility','off');
plot(ax,[x x],[yLow yHigh],'-','Color',edge,'LineWidth',1.05,'HandleVisibility','off');
plot(ax,x,y,marker,'LineStyle','none','MarkerSize',size,'MarkerFaceColor',color, ...
    'MarkerEdgeColor',edge,'LineWidth',.8,'HandleVisibility','off');
end

function rgb=hex2rgb(hex)
hex=erase(string(hex),'#');
rgb=[hex2dec(extractBetween(hex,1,2)),hex2dec(extractBetween(hex,3,4)), ...
    hex2dec(extractBetween(hex,5,6))]/255;
end
