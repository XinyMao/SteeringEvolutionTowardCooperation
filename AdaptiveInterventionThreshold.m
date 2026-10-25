function [fig, summaryTable, landscapeData] = plot_fig5c_final(localCsv, pcCsv, globalCsv, saveBase)
%PLOT_FIG5C_FINAL Nature-style four-axis Fig. 5c.
% Three local-drift landscapes (beta=0,2,5) plus threshold summary.
% This function only reads existing simulation outputs; it does not
% re-estimate thresholds or alter analytical formulas.
if nargin < 1 || isempty(localCsv)
    localCsv=fullfile('results','fig5c_beta_x_supplement','fig5c_local_drift.csv');
end
if nargin < 2 || isempty(pcCsv)
    pcCsv=fullfile('results','fig5c_beta_x_supplement','fig5c_pc_by_x.csv');
end
if nargin < 3 || isempty(globalCsv)
    globalCsv=fullfile('results','fig5c_beta_x_supplement','fig5c_global_threshold.csv');
end
if nargin < 4, saveBase=''; end

L=readtable(localCsv,'TextType','string');
P=readtable(pcCsv,'TextType','string');
G=readtable(globalCsv,'TextType','string');
betas=[0 2 5]; k=10;
xv=unique(L.x_target,'sorted'); pv=unique(L.p_target,'sorted');

% Aggregate mean drift on the original x/p grid for the three requested beta values.
M=cell(numel(betas),1); theoryLocal=cell(numel(betas),1);
for ib=1:numel(betas)
    D=L(abs(L.beta-betas(ib))<1e-12,:); Z=nan(numel(pv),numel(xv));
    for ip=1:numel(pv)
        for ix=1:numel(xv)
            q=D(abs(D.p_target-pv(ip))<1e-12 & abs(D.x_target-xv(ix))<1e-12,:);
            if ~isempty(q), Z(ip,ix)=mean(q.drift,'omitnan'); end
        end
    end
    M{ib}=Z;
    q=P(abs(P.beta-betas(ib))<1e-12,:); q=sortrows(q,'x_target');
    theoryLocal{ib}=q;
end
allVals=[]; for ib=1:numel(M), allVals=[allVals; M{ib}(isfinite(M{ib}))]; end %#ok<AGROW>
Dmax=max(abs(allVals));

fig=figure('Color','w','Name','Fig. 5c final','Units','normalized','Position',[.04 .12 .92 .72]);
% Four ordinary axes: all positions remain directly editable in MATLAB.
leftPos={[.055 .20 .155 .66],[.235 .20 .155 .66],[.415 .20 .155 .66]};
axs=gobjects(1,4);
cmap=driftMap(256); colormap(fig,cmap);
for ib=1:numel(betas)
    ax=axes('Parent',fig,'Position',leftPos{ib}); axs(ib)=ax;
    % Smooth visual rendering of the mean-drift grid only; all boundaries
    % and theory curves remain based on the original simulation grid/data.
    xq=linspace(min(xv),max(xv),220); pq=linspace(min(pv),max(pv),220);
    try, Mq=interp2(xv,pv,M{ib},xq,pq','makima');
    catch, Mq=interp2(xv,pv,M{ib},xq,pq','spline'); end
    imagesc(ax,xq,pq,Mq); set(ax,'YDir','normal'); caxis(ax,[-Dmax Dmax]);
    xlim(ax,[min(xv) max(xv)]); ylim(ax,[min(pv) max(pv)]); hold(ax,'on'); box(ax,'off');
    % Empirical zero-drift contour (black solid).
    contour(ax,xq,pq,Mq,[0 0],'k-','LineWidth',1.25);
    % Theoretical local threshold (white dashed), without any displacement.
    q=theoryLocal{ib}; plot(ax,q.x_target,q.p_c_theory,'w--','LineWidth',1.35);
    text(ax,min(xv)+.035*range(xv),max(pv)-.075*range(pv),sprintf('\\beta = %g',betas(ib)), ...
        'Color','k','FontName','Arial','FontSize',10,'FontWeight','bold','Interpreter','tex');
    xlabel(ax,'Initial cooperation x');
    if ib==1, ylabel(ax,'Steering fraction p'); else, set(ax,'YTickLabel',[]); end
    set(ax,'FontName','Arial','FontSize',9,'Layer','top','TickDir','out');
end
% One shared colorbar for the three landscapes.
cb=colorbar('peer',axs(3)); cb.Position=[.585 .20 .014 .66]; cb.Label.String='Local evolutionary drift';

% Threshold summary: retain all available beta values in the right panel.
axs(4)=axes('Parent',fig,'Position',[.665 .20 .285 .66]); ax=axs(4); hold(ax,'on'); box(ax,'off');
G=sortrows(G,'beta');
% finite-grid theory envelope over exactly the simulated x grid
finiteTheory=nan(height(G),1); controller=nan(height(G),1);
for ib=1:height(G)
    q=P(abs(P.beta-G.beta(ib))<1e-12,:);
    finiteTheory(ib)=max(q.p_c_theory);
    [~,jj]=max(q.p_c_empirical);
    controller(ib)=q.x_target(jj);
end
pstatic=G.p_static_star(1);
finiteTheory=finiteTheory./pstatic;
% Draw low-priority theory first.
plot(ax,G.beta,finiteTheory,'--','Color',[.45 .45 .45],'LineWidth',1.05, ...
    'DisplayName','Theory, simulated x-grid');
plot(ax,G.beta,exp(-G.beta./k),'-','Color',[0.02 0.25 0.52],'LineWidth',1.8, ...
    'DisplayName','Global theory');
errorbar(ax,G.beta,G.normalized_empirical, ...
    G.normalized_empirical-G.empirical_ci95_low./pstatic, ...
    G.empirical_ci95_high./pstatic-G.normalized_empirical,'o', ...
    'Color',[.05 .05 .05],'MarkerFaceColor',[.05 .05 .05],'MarkerSize',5.5, ...
    'LineWidth',1.15,'DisplayName','Simulation envelope');
xlabel(ax,'Targeting strength \\beta','Interpreter','tex');
ylabel(ax,'Normalized steering threshold p / p_{static}^*','Interpreter','tex');
title(ax,'Threshold reduction','FontWeight','normal');
legend(ax,'Location','southwest','Box','off'); grid(ax,'on'); ax.GridAlpha=.10;
set(ax,'FontName','Arial','FontSize',9,'Layer','top','TickDir','out');
annotation(fig,'textbox',[.018 .91 .03 .05],'String','c','EdgeColor','none', ...
    'FontName','Arial','FontSize',14,'FontWeight','bold','HorizontalAlignment','center');

summaryTable=table(G.beta,G.normalized_empirical,G.empirical_ci95_low./pstatic, ...
    G.empirical_ci95_high./pstatic,finiteTheory,exp(-G.beta./k),controller, ...
    'VariableNames',{'beta','empirical_envelope','empirical_ci95_low', ...
    'empirical_ci95_high','finite_grid_theory','global_analytical_theory','controlling_x'});
landscapeData=struct('beta',betas,'x',xv,'p',pv,'meanDrift',{M}, ...
    'theoryLocal',{theoryLocal},'Dmax',Dmax,'axes',axs,'colorbar',cb);
if ~isempty(saveBase)
    savefig(fig,[saveBase,'.fig']);
    exportgraphics(fig,[saveBase,'.pdf'],'ContentType','vector');
    writetable(summaryTable,[saveBase,'_summary.csv']);
    try, exportgraphics(fig,[saveBase,'.svg'],'ContentType','vector'); catch, end
end
end

function cmap=driftMap(n)
% Requested muted diverging palette: negative blue, zero near-white, positive red.
anchors=[0.0235 0.4431 0.6902; 0.5725 0.7725 0.8706; 0.9686 0.9686 0.9686; ...
         0.9569 0.6471 0.5098; 0.7922 0.0000 0.1255];
t=linspace(0,1,n); ta=linspace(0,1,size(anchors,1));
cmap=[interp1(ta,anchors(:,1),t)' interp1(ta,anchors(:,2),t)' interp1(ta,anchors(:,3),t)'];
end
