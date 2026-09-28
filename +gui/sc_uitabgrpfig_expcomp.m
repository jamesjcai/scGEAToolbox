function sc_uitabgrpfig_expcomp(sce, glist, parentfig, ~, thisc)

if nargin < 3, parentfig = []; end

% https://www.mathworks.com/help/rptgen/ug/compile-a-presentation-program.html
if (ismcc || isdeployed) && pkg.i_isreportgenavailable('ppt'), makePPTCompilable(); end


hx=gui.myFigure(parentfig, true);
hFig=hx.FigHandle;
hFig.Position(3) = hFig.Position(3) * 1.8;

% if ~isempty(parentfig) && isa(parentfig,'matlab.ui.Figure')
%     px = hFig.Position;
%     px_new = [cx(1)-px(3)/2 cx(2)-px(4)/2];
%
%     % if px_new(1)<0
%     %     ss = get(0, 'screensize');
%     %     px_new(1)=px_new(1)-ss(4); %ss(3);
%     % end
% else
%     px_new = [];
% end

n = length(glist);

tabgp = uitabgroup();
tab = cell(n,1);
ax0 = cell(n,1);
ax = cell(n,2);

idx = 1;
focalg = glist(idx);

[c, cL] = findgroups(string(thisc));
sce1=copy(sce);
sce1 = sce1.selectcells(c == 1);  % OK
sce2=copy(sce);
sce2 = sce2.selectcells(c == 2);  % OK

for k=1:n
    c1 = sce1.X(sce1.g == glist(k), :);
    if issparse(c1), c1 = full(c1); end
    c2 = sce2.X(sce2.g == glist(k), :);
    if issparse(c2), c2 = full(c2); end

    tab{k} = uitab(tabgp, 'Title', sprintf('%s',glist(k)));

    ax0{k} = axes('parent',tab{k});
    ax{k,1} = subplot(1,2,1);

    scatter(sce1.s(:,1), sce1.s(:,2), 5, c1, 'filled');
    stem3(sce1.s(:,1), sce1.s(:,2), c1, 'marker', 'none', 'color', 'm');
    hold on;
    scatter3(sce1.s(:,1), sce1.s(:,2), zeros(size(sce1.s(:,2))), 5, c1, 'filled');
    z1 = zlim(ax{k,1});
    z1 = z1(2);


    ax{k,2} = subplot(1,2,2);
    scatter(sce2.s(:,1), sce2.s(:,2), 5, c2, 'filled');
    stem3(sce2.s(:,1), sce2.s(:,2), c2, 'marker', 'none', 'color', 'm');
    hold on;
    scatter3(sce2.s(:,1), sce2.s(:,2), zeros(size(sce2.s(:,2))), 5, c2, 'filled');
    z2 = zlim(ax{k,2});
    z2 = z2(2);

    zz = max([z1 z2]);
    zlim(ax{k,1},[0 zz]);
    zlim(ax{k,2},[0 zz]);


    title(ax{k,1}, glist(k)+" "+string(cL{1}));
    subtitle(ax{k,1}, gui.i_getsubtitle(c1));
    title(ax{k,2}, glist(k)+" "+string(cL{2}));
    subtitle(ax{k,2}, gui.i_getsubtitle(c2));
end

tabgp.SelectionChangedFcn=@displaySelection;

hx.addCustomButton('off',  @i_genecards, 'www.jpg', 'GeneCards...');
hx.addCustomButton('off', @i_proteinstructure, 'hexagon_16dp_000000_FILL0_wght400_GRAD0_opsz20.jpg', 'Protein Structure...');
hx.addCustomButton('off', @in_savedata, "Save.gif", 'Save Gene List...');

hx.show(parentfig);


function in_savedata(~,~)
        gui.i_exporttable(table(glist), true, ...
            'Tmarkerlist', 'MarkerListTable', [], [], hFig);
    end

function displaySelection(~,event)
        t = event.NewValue;
        txt = t.Title;
        [~,idx]=ismember(txt,glist);
        focalg = glist(idx);
    end

function i_genecards(~, ~)
        web(sprintf('https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s', focalg),'-new');
    end

function i_proteinstructure(~, ~)
        gui.i_viewprotein(focalg, ParentFig=hFig);
    end

end
