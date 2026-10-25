function [fig, phaseData, snapshotInfo] = plot_fig5d_final(phaseCsv, gamma0Dir, gamma20Dir, saveBase, replicate, times)
%PLOT_FIG5D_FINAL Nature Communications-style Fig. 5d.
if nargin < 1 || isempty(phaseCsv), phaseCsv=fullfile('results','fig5d_phase_space_formal_v3','fig5d_trajectories.csv'); end
if nargin < 2 || isempty(gamma0Dir), gamma0Dir=fullfile('results','fig5d_animation_x040_gamma0_r5','snapshots_mat'); end
if nargin < 3 || isempty(gamma20Dir), gamma20Dir=fullfile('results','fig5d_animation_x040_gamma20_r5','snapshots_mat'); end
if nargin < 4, saveBase=''; end
if nargin < 5, replicate=0; end
if nargin < 6 || isempty(times), times=[0 40]; end
if numel(times)~=2 || times(1)~=0 || times(2)~=40, error('Use exactly times [0 40].'); end

T=readtable(phaseCsv,'TextType','string');
req={'x','z','time_mcs','beta','p','k','z_ad','z_rand','z_n_star','replicate'};
if any(~ismember(req,T.Properties.VariableNames)), error('Missing phase-space columns.'); end
if ~ismember('x0',T.Properties.VariableNames), T.x0=recoverX0(T); end
p=mean(T.p); k=mean(T.k); D=T(abs(T.beta-.5)<1e-12 & abs(T.p-.03)<1e-12,:); x0v=unique(D.x0,'sorted');
fig=figure('Color','w','Name','Fig. 5d final','Units','normalized','Position',[.03 .08 .94 .82]);
axPhase=axes('Parent',fig,'Position',[.055 .22 .36 .64]); hold(axPhase,'on'); box(axPhase,'off');
xg=linspace(0,1,600)'; a=p+(1-p).*xg; zad=a-2*a.^2.*(1-a)./(a.*(1-p)+p); zrand=a.^2;
den=(k-1).*(a.*(1+p)-p); zn=nan(size(xg)); ok=den>0; zn(ok)=a(ok).^2.*(1+(k-2).*a(ok))./den(ok);
patch(axPhase,[xg;flipud(xg)],[zad;ones(size(zad))],[.92 .95 .96],'FaceAlpha',.10,'EdgeColor','none','HandleVisibility','off');
h1=plot(axPhase,xg,zad,'-','Color',[.08 .10 .12],'LineWidth',2.0,'DisplayName','z_{ad}(x)');
h2=plot(axPhase,xg,zn,'-.','Color',[.27 .38 .40],'LineWidth',1.2,'DisplayName','z_n^*(x)');
h3=plot(axPhase,xg,zrand,':','Color',[.62 .64 .64],'LineWidth',1.0,'DisplayName','z_{rand}(x)');
cols=[.32 .43 .51; .38 .53 .56; .49 .60 .55; .65 .56 .43; .67 .43 .34]; phaseData=struct();
for ix=1:numel(x0v)
    Q=D(abs(D.x0-x0v(ix))<1e-12,:); reps=unique(Q.replicate,'sorted'); light=.45+.55*cols(ix,:); mt=unique(Q.time_mcs,'sorted'); mx=nan(size(mt)); mz=mx;
    for ir=1:numel(reps)
        R=sortrows(Q(Q.replicate==reps(ir),:),'time_mcs'); plot(axPhase,R.x,R.z,'-','Color',light,'LineWidth',.65,'HandleVisibility','off');
    end
    for it=1:numel(mt), R=Q(abs(Q.time_mcs-mt(it))<1e-12,:); mx(it)=mean(R.x); mz(it)=mean(R.z); end
    plot(axPhase,mx,mz,'-','Color',cols(ix,:),'LineWidth',1.05,'HandleVisibility','off');
    plot(axPhase,mx(1),mz(1),'o','Color',cols(ix,:),'MarkerFaceColor','w','MarkerSize',5,'HandleVisibility','off');
    plot(axPhase,mx(end),mz(end),'o','Color',cols(ix,:),'MarkerFaceColor',cols(ix,:),'MarkerSize',5,'HandleVisibility','off');
    text(axPhase,mx(1)+.012,mz(1),sprintf('x_0=%.1f',x0v(ix)),'Color',cols(ix,:),'FontSize',7,'FontName','Arial','HandleVisibility','off');
    phaseData(ix).x0=x0v(ix); phaseData(ix).time=mt; phaseData(ix).x_mean=mx; phaseData(ix).z_mean=mz;
end
text(axPhase,.60,.84,'Targeting-beneficial','Color',[.33 .40 .42],'FontSize',8,'FontName','Arial');
xlabel(axPhase,'Cooperation level x'); ylabel(axPhase,'Pair density z = p_{AA}'); xlim(axPhase,[0 1]); ylim(axPhase,[0 1]);
legend(axPhase,[h1 h2 h3],'Location','southeast','Box','off'); set(axPhase,'FontName','Arial','FontSize',10,'Layer','top','TickDir','out');

dInitial=pickSnap(gamma20Dir,.5,replicate,0); G=graph(dInitial.edges(:,1),dInitial.edges(:,2));
try, tmp=figure('Visible','off'); pp=plot(G,'Layout','force','NodeLabel',{}); xy=[pp.XData(:),pp.YData(:)]; close(tmp);
catch, th=linspace(0,2*pi,numnodes(G)+1)'; xy=[cos(th(1:end-1)),sin(th(1:end-1))]; end
dFixed=pickSnap(gamma0Dir,0,replicate,40); dUnbiased=pickSnap(gamma20Dir,0,replicate,40); dTargeted=pickSnap(gamma20Dir,.5,replicate,40);
snapshots={dInitial,dFixed,dUnbiased,dTargeted}; labels={'Initial','Fixed','Unbiased','Targeted'};
nbh=neighbors(G,1); zoomNodes=unique([1; nbh(:)]); bnd=zoomBounds(xy,zoomNodes); sx=[.475 .590 .705 .820]; axs=gobjects(1,4);
for j=1:4
    axs(j)=axes('Parent',fig,'Position',[sx(j) .56 .105 .30]); drawGlobal(axs(j),G,xy,snapshots{j},j>=3,bnd);
    title(axs(j),{labels{j},sprintf('t = %.0f',snapshots{j}.time_mcs)},'FontSize',8,'FontWeight','normal');
end
zoomAxes=gobjects(1,2);
for j=1:2
    zoomAxes(j)=axes('Parent',fig,'Position',[sx(j+2) .19 .105 .25]); drawZoom(zoomAxes(j),G,xy,snapshots{j+2},bnd); title(zoomAxes(j),labels{j+2},'FontSize',8,'FontWeight','normal');
end
annotation(fig,'textbox',[.018 .91 .03 .05],'String','d','EdgeColor','none','FontName','Arial','FontSize',14,'FontWeight','bold','HorizontalAlignment','center');

sameGraph=isequal(dInitial.edges,dFixed.edges)&&isequal(dInitial.edges,dUnbiased.edges)&&isequal(dInitial.edges,dTargeted.edges);
sameStrategy=isequal(dInitial.strategy,dFixed.strategy)&&isequal(dInitial.strategy,dUnbiased.strategy)&&isequal(dInitial.strategy,dTargeted.strategy);
sameSteering=isequal(dInitial.is_steering,dFixed.is_steering)&&isequal(dInitial.is_steering,dUnbiased.is_steering)&&isequal(dInitial.is_steering,dTargeted.is_steering);
snapshotInfo=struct('N',double(dInitial.N),'k',double(dInitial.k),'p',double(dInitial.p),'x0',.4,'replicate',replicate,'times',times,'networkSeed',double(dInitial.seed), ...
    'sameInitialGraph',sameGraph,'sameInitialStrategy',sameStrategy,'sameInitialSteering',sameSteering,'fixedNodeCoordinates',true,'zoomFocalNode',1,'zoomNodeIDs',zoomNodes,'zoomBounds',bnd, ...
    'phaseBeta',.5,'phaseP',.03,'phaseAxes',axPhase,'snapshotAxes',axs,'zoomAxes',zoomAxes,'nodeCoordinates',xy,'cooperatorRGB',[.36 .53 .64],'defectorRGB',[.72 .47 .34], ...
    'steeringRingRGB',[.34 .70 .28],'edgeRGB',[.82 .82 .82],'strategyMarkerSize',9,'steeringRingSize',18,'steeringRingLineWidth',1.35);
if ~isempty(saveBase), savefig(fig,[saveBase,'.fig']); exportgraphics(fig,[saveBase,'.pdf'],'ContentType','vector'); try, exportgraphics(fig,[saveBase,'.svg'],'ContentType','vector'); catch, end, end
end

function x0=recoverX0(T)
x0=nan(height(T),1); starts=find(abs(T.time_mcs-min(T.time_mcs))<1e-12);
for j=1:numel(starts), a=starts(j); if j<numel(starts), b=starts(j+1)-1; else, b=height(T); end; x0(a:b)=T.x(a); end
if any(~isfinite(x0)), error('Could not recover x0.'); end
end
function d=pickSnap(dirName,beta,rep,tWanted)
ff=dir(fullfile(dirName,sprintf('*x0.40000_b%.3f_rep%04d*.mat',beta,rep))); if isempty(ff), error('No RRG snapshot for beta=%.3g, replicate=%d.',beta,rep); end
tt=zeros(numel(ff),1); for j=1:numel(ff), q=load(fullfile(ff(j).folder,ff(j).name),'time_mcs'); tt(j)=double(q.time_mcs); end
[~,j]=min(abs(tt-tWanted)); d=load(fullfile(ff(j).folder,ff(j).name));
end
function drawGlobal(ax,G,xy,d,mark,bnd)
plot(ax,G,'XData',xy(:,1),'YData',xy(:,2),'NodeLabel',{},'Marker','none','EdgeColor',[.84 .84 .84],'LineWidth',.10); hold(ax,'on'); drawNodes(ax,xy,d,9,18,1.35);
if mark, rectangle(ax,'Position',[bnd(1) bnd(2) bnd(3)-bnd(1) bnd(4)-bnd(2)],'EdgeColor',[.20 .20 .20],'LineWidth',.7); end; axis(ax,'equal'); axis(ax,'off');
text(ax,.02,.03,sprintf('x = %.3f, z = %.3f',d.x,d.z),'Units','normalized','FontName','Arial','FontSize',7,'Color',[.20 .20 .20]);
end
function drawZoom(ax,G,xy,d,bnd)
plot(ax,G,'XData',xy(:,1),'YData',xy(:,2),'NodeLabel',{},'Marker','none','EdgeColor',[.72 .72 .72],'LineWidth',.22); hold(ax,'on'); drawNodes(ax,xy,d,12,22,1.55); axis(ax,'equal'); axis(ax,'off');
pad=.12*max(bnd(3)-bnd(1),bnd(4)-bnd(2)); xlim(ax,[bnd(1)-pad bnd(3)+pad]); ylim(ax,[bnd(2)-pad bnd(4)+pad]);
text(ax,.03,.04,sprintf('x = %.3f, z = %.3f',d.x,d.z),'Units','normalized','FontName','Arial','FontSize',7,'Color',[.20 .20 .20]);
end
function drawNodes(ax,xy,d,nodeSize,ringSize,ringWidth)
s=logical(d.is_steering(:)); c=~s & logical(d.strategy(:)); q=~s & ~logical(d.strategy(:));
scatter(ax,xy(q,1),xy(q,2),nodeSize,[.72 .47 .34],'filled','MarkerFaceAlpha',.58); scatter(ax,xy(c,1),xy(c,2),nodeSize,[.36 .53 .64],'filled','MarkerFaceAlpha',.58);
scatter(ax,xy(s,1),xy(s,2),ringSize,'o','MarkerFaceColor','none','MarkerEdgeColor','w','LineWidth',1.8); scatter(ax,xy(s,1),xy(s,2),ringSize,'o','MarkerFaceColor','none','MarkerEdgeColor',[.34 .70 .28],'LineWidth',ringWidth);
end
function b=zoomBounds(xy,nodes)
q=xy(nodes,:); r=max(max(q,[],1)-min(q,[],1)); if r==0, r=.05; end; m=.18*r; b=[min(q(:,1))-m,min(q(:,2))-m,max(q(:,1))+m,max(q(:,2))+m];
end
