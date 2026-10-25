function [fig, gridData] = plot_fig5b_exposure_heatmap(dataPath, saveBase)
%PLOT_FIG5B_EXPOSURE_HEATMAP Exposure landscape and beta slices.
if nargin < 1 || isempty(dataPath)
    dataPath=fullfile('results','fig5b_formal_exposure','fig5b_summary.csv');
end
if nargin < 2, saveBase=''; end
T=readtable(dataPath,'TextType','string');
xv=unique(T.x_target,'sorted'); bv=unique(T.beta,'sorted');
Z=nan(numel(xv),numel(bv)); Zlo=Z; Zhi=Z; Th=Z;
for ix=1:numel(xv)
    for ib=1:numel(bv)
        use=abs(T.x_target-xv(ix))<1e-12 & abs(T.beta-bv(ib))<1e-12;
        y=T.eta_normalized(use); y=y(isfinite(y));
        if isempty(y), continue; end
        Z(ix,ib)=mean(y); se=std(y,0)/sqrt(numel(y));
        Zlo(ix,ib)=Z(ix,ib)-1.96*se; Zhi(ix,ib)=Z(ix,ib)+1.96*se;
        q=T.theory_eta_normalized(use); q=q(isfinite(q));
        if ~isempty(q), Th(ix,ib)=mean(q); end
    end
end
fig=figure('Color','w','Name','Fig. 5b exposure landscape');
% Use ordinary axes rather than tiledlayout so panels remain freely movable
% and resizable during final multi-panel composition.
ax=axes('Parent',fig,'Position',[0.08 0.16 0.38 0.72]);
G=100*(Z-1); % percent gain relative to beta=0
imagesc(ax,1:numel(bv),1:numel(xv),G); set(ax,'YDir','normal'); axis(ax,'tight');
gainMap=blueTealOrange(256);
colormap(ax,gainMap); clim(ax,[0 max(G(:),[],'omitnan')]);
cb=colorbar(ax); cb.Label.String='Relative exposure gain (%)';
xlabel(ax,'Targeting strength \beta'); ylabel(ax,'Initial cooperation x');
set(ax,'XTick',1:numel(bv),'XTickLabel',compose('%.2g',bv), ...
    'YTick',1:numel(xv),'YTickLabel',compose('%.2g',xv));
title(ax,'Relative exposure gain'); hold(ax,'on');
for ix=1:numel(xv)
    for ib=1:numel(bv)
        if isfinite(G(ix,ib))
            u=(G(ix,ib)-0)/(max(G(:),[],'omitnan')-0+eps);
            ci=max(1,min(size(gainMap,1),1+round(u*(size(gainMap,1)-1))));
            rgb=gainMap(ci,:); lum=0.299*rgb(1)+0.587*rgb(2)+0.114*rgb(3);
            % Article display rule: gains up to and including 10% use white;
            % larger gains use black for the lighter warm background.
            txtColor='w'; if G(ix,ib)>10, txtColor='k'; end
            text(ax,ib,ix,sprintf('+%.1f%%',G(ix,ib)), ...
                'HorizontalAlignment','center','FontSize',8,'Color',txtColor);
        end
    end
end
grid(ax,'on'); ax.GridAlpha=.15;
ax2=axes('Parent',fig,'Position',[0.56 0.16 0.38 0.72]);
hold(ax2,'on'); box(ax2,'on'); cols=lines(numel(xv));
for ix=1:numel(xv)
    errorbar(ax2,bv,Z(ix,:),Z(ix,:)-Zlo(ix,:),Zhi(ix,:)-Z(ix,:),'o-', ...
        'Color',cols(ix,:),'MarkerFaceColor',cols(ix,:),'LineWidth',1.4, ...
        'DisplayName',sprintf('x = %.2g',xv(ix)));
    % Draw theory after simulation so a close curve remains visible.
    if ix==1, theoryName='weak-targeting theory'; theoryVisibility='on';
    else, theoryName=''; theoryVisibility='off'; end
    plot(ax2,bv,Th(ix,:),'--','Color',cols(ix,:),'LineWidth',1.6, ...
        'DisplayName',theoryName,'HandleVisibility',theoryVisibility);
end
yline(ax2,1,':','Color',[.3 .3 .3],'HandleVisibility','off');
xlabel(ax2,'Targeting strength \beta'); ylabel(ax2,'Normalized exposure');
title(ax2,'Slices with 95% CI and theory'); legend(ax2,'Location','northwest','Box','off');
grid(ax2,'on'); ax2.GridAlpha=.15; set([ax ax2],'FontName','Arial','FontSize',10,'Layer','top');
tx=repmat(xv(:),numel(bv),1);
tb=repmat(bv(:),numel(xv),1);
tt=reshape(Th.',[],1);
theoryTable=table(tx,tb,tt,'VariableNames',{'x','beta','theory_eta_normalized'});
gridData=struct('x',xv,'beta',bv,'mean',Z,'ci95_low',Zlo,'ci95_high',Zhi,'theory',Th, ...
    'theoryTable',theoryTable,'axHeatmap',ax,'axCurves',ax2);
if ~isempty(saveBase), savefig(fig,[saveBase,'.fig']); exportgraphics(fig,[saveBase,'.pdf'],'ContentType','vector'); end
end

function cmap=blueTealOrange(n)
% Palette sampled from color.png: dark blue/slate -> orange -> muted teal.
% The mid-range orange and high-range teal are intentional features of the
% article's existing color axis; do not replace this with jet/rainbow.
anchors = [0.28 0.40 0.52
            0.38 0.55 0.63   % blue-teal
            0.52 0.70 0.68   % muted cyan
            0.66 0.74 0.64   % sage
            0.80 0.72 0.52   % muted gold
            0.84 0.58 0.45   % soft terracotta
        ];
t=linspace(0,1,n); ta=linspace(0,1,size(anchors,1));
cmap=[interp1(ta,anchors(:,1),t)' interp1(ta,anchors(:,2),t)' ...
       interp1(ta,anchors(:,3),t)'];
end
