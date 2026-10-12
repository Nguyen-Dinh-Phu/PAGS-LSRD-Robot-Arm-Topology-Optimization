
%----------------------------- PAGS_Output.m -----------------------------%
% Input : mode    - Output mode ('history' or 'sweep')
%         results - Optimization results
%         mesh    - Finite element mesh and geometry
%         outDir  - Output directory
%-------------------------------------------------------------------------%

function PAGS_Output(mode,results,mesh,~)
switch lower(char(mode))
    case 'history'
        baselineEvolution(results.baseline,mesh);             
        baselineTopology(results.baseline,mesh);             
        convergencePlots(results.baseline,results.angles);    
    case 'sweep'
        topologyMatrix(results,mesh);                         
        tradeoffPlots(results);                             
end
end

%----------------------------- BASELINE EVOLUTION
function baselineEvolution(base,mesh)
fig = newFigure('Fig. 3 - Baseline evolution',[90 100 1200 600]);
ax = axes('Parent',fig,'Position',[.09 .14 .82 .74]);
yyaxis(ax,'left');
plot(ax,base.h_iter,base.h_vol,'-ks','LineWidth',1.3, ...
    'MarkerSize',3.5,'MarkerIndices',1:10:numel(base.h_iter));
ylabel(ax,'Volume fraction');
yyaxis(ax,'right');
plot(ax,base.h_iter,base.h_comp,'-ro','LineWidth',1.3, ...
    'MarkerSize',3.5,'MarkerIndices',1:10:numel(base.h_iter));
ylabel(ax,'Compliance');
xlabel(ax,'Iteration'); xlim(ax,[0 400]);
grid(ax,'on'); box(ax,'on');
set(ax,'FontName','Times New Roman','FontSize',12);
legend(ax,{'Volume fraction','Compliance'},'Location','northwest');
pos = [.11 .62 .15 .16; .28 .53 .15 .16; .44 .43 .15 .16; ...
       .61 .33 .15 .16; .77 .22 .15 .16];
for k = 1:numel(base.snapshotIter)
    a = axes('Parent',fig,'Position',pos(k,:));
    drawMaterial(a,base.snapshotPhi{k},mesh);
    title(a,sprintf('k=%d',base.snapshotIter(k)), ...
        'FontName','Times New Roman','FontSize',9);
end
showFigure(fig);
end

%----------------------------- BASELINE TOPOLOGY
function baselineTopology(base,mesh)
fig = newFigure('Fig. 4 - Baseline topology',[100 100 1100 390]);
ax = axes('Parent',fig,'Position',[.05 .10 .91 .80]);
hold(ax,'on');
phi = reshape(base.phi_final,mesh.nely+1,mesh.nelx+1);
contourf(ax,mesh.X,mesh.Y,phi,[0 max(phi(:))+1], ...
    'FaceColor',[.15 .30 .45],'LineColor','none');
contour(ax,mesh.X,mesh.Y,mesh.Y, ...
    min(mesh.Y(:)):1:max(mesh.Y(:)), ...
    'Color',[.85 .3 .2],'LineWidth',.65);        % Fiber direction guides
contourf(ax,mesh.X,mesh.Y,-phi,[0 max(-phi(:))+1], ...
    'FaceColor','w','LineColor','none');         % Void regions
contour(ax,mesh.X,mesh.Y,phi,[0 0],'k','LineWidth',.5);
axis(ax,'equal','tight','off');
showFigure(fig);
end

%----------------------------- CONVERGENCE PLOTS
function convergencePlots(base,angles)
names = {'Fig. 5a - Volume fraction','Fig. 5b - Compliance', ...
         'Fig. 6a - CG-RC distance','Fig. 6b - Objective'};
fields = {'h_vol','h_comp','h_d2','h_obj'};
labels = {'Volume fraction','Compliance', ...
          'Squared CG-RC distance (mm^2)','Objective function J'};
colors = {'k','b','g','m','r'};
markers = {'*','o','^','s','x'};
cases = [base,angles];
for m = 1:4
    fig = newFigure(names{m},[100 100 930 650]);
    ax = axes('Parent',fig);
    setupAxes(ax,'Iteration',labels{m});
    for j = 1:numel(cases)
        c = cases(j);
        if j == 1
            name = 'Base (w=0, \theta=0^\circ)';
        else
            name = sprintf('w=%.1f, \\theta=%g^\\circ',c.w,c.theta_deg);
        end
        plot(ax,c.h_iter,c.(fields{m}),'-', ...
            'Color',colors{j},'Marker',markers{j}, ...
            'LineWidth',1.6,'MarkerSize',5, ...
            'MarkerIndices',1:10:numel(c.h_iter), ...
            'DisplayName',name);
    end

    xlim(ax,[0 400]);
    legend(ax,'show','Location','best','FontSize',11);
    showFigure(fig);
end
end
%----------------------------- TOPOLOGY MATRIX
function topologyMatrix(results,mesh)
[nAngles,nWeights] = size(results);
fig = newFigure('Table 3 - Optimized topologies',[80 50 1450 980]);
tl = tiledlayout(fig,nWeights,nAngles, ...
    'TileSpacing','compact','Padding','loose');
for iw = 1:nWeights
    for ia = 1:nAngles
        c = results(ia,iw);
        ax = nexttile(tl,(iw-1)*nAngles+ia);
        drawMaterial(ax,c.phi_final,mesh);
        if iw == 1
            title(ax,sprintf('\\theta=%g^\\circ',c.theta_deg), ...
                'FontName','Times New Roman','FontSize',15);
        end
        if ia == 1
            text(ax,-.2*mesh.fem.domain.lx,mesh.fem.domain.ly/2, ...
                sprintf('w=%.1f',c.w), ...
                'FontName','Times New Roman','FontSize',13, ...
                'Clipping','off','HorizontalAlignment','right');
        end
    end
end
showFigure(fig);
end

%----------------------------- TRADE-OFF PLOTS
function tradeoffPlots(results)
names = {'Fig. 7a - Normalized distance', ...
         'Fig. 7b - Normalized compliance','Fig. 7c - Trade-off'};
xLabels = {'PAGS weight, w','PAGS weight, w','C/C_0'};
yLabels = {'D/D_0','C/C_0','D/D_0'};
colors = {'b','g','m','r'};
markers = {'o','^','s','x'};
[nAngles,nWeights] = size(results);
for m = 1:3
    fig = newFigure(names{m},[80 80 850 630]);
    ax = axes('Parent',fig);
    setupAxes(ax,xLabels{m},yLabels{m});
    for ia = 1:nAngles
        c = results(ia,:);
        w = [c.w]; C = [c.Cnorm_final]; D = [c.Dnorm_final];
        switch m
            case 1, x = w; y = D;
            case 2, x = w; y = C;
            case 3, x = C; y = D;
        end
        plot(ax,x,y,'-','Color',colors{ia}, ...
            'Marker',markers{ia},'LineWidth',1.6,'MarkerSize',6, ...
            'DisplayName',sprintf('\\theta=%g^\\circ',c(1).theta_deg));
        if m == 3
            for iw = 1:nWeights
                text(ax,x(iw),y(iw),sprintf(' %.1f',w(iw)), ...
                    'FontName','Times New Roman','FontSize',9);
            end
        end
    end
    if m < 3, xlim(ax,[0 .5]); end
    legend(ax,'show','Location','best','FontSize',11);
    showFigure(fig);
end
end

%----------------------------- PLOTTING UTILITIES
function drawMaterial(ax,phiNodal,mesh)
phi = reshape(phiNodal,mesh.nely+1,mesh.nelx+1);
contourf(ax,mesh.X,mesh.Y,phi,[0 max(phi(:))+1], ...
    'FaceColor','k','LineColor','none');
axis(ax,'equal','tight','off');
end

%-----------------------------
function setupAxes(ax,xLabel,yLabel)
grid(ax,'on'); box(ax,'on'); hold(ax,'on');
set(ax,'FontName','Times New Roman','FontSize',13,'LineWidth',.8);
xlabel(ax,xLabel,'Interpreter','tex');
ylabel(ax,yLabel,'Interpreter','tex');
end

%-----------------------------
function fig = newFigure(name,pos)
fig = figure('Name',name,'NumberTitle','off', ...
    'Color','w','Visible','off','Position',pos);
end

%-----------------------------
function showFigure(fig)
set(fig,'Visible','on');                       
end
%-------------------------------------------------------------------------%
