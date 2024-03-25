function callback_DVGene2Groups(src, ~)

%warndlg('The function is under development.');
%return;

    FigureHandle = src.Parent.Parent;
    sce = guidata(FigureHandle);
    if ~gui.gui_showrefinfo('DV Analysis'), return; end
    
    [i1, i2, cL1, cL2] = gui.i_select2grps(sce);
    if length(i1) == 1 || length(i2) == 1, return; end
    
    c=zeros(size(i1));
    c(i1)=1; c(i2)=2;
    cL=[cL1;cL2];
    cL1 = matlab.lang.makeValidName(cL1);
    cL2 = matlab.lang.makeValidName(cL2);
    if ~all(c>0)
        sce=sce.selectcells(c>0);
        c=c(c>0);
        i1=c==1;
        i2=c==2;
    end
    
    sce1 = sce.selectcells(i1);
    sce2 = sce.selectcells(i2);
    sce1 = sce1.qcfilter;  
    sce2 = sce2.qcfilter;  

    if sce1.NumCells < 50 || sce2.NumCells < 50
        uiwait(warndlg('One of groups contains too few cells (n < 50). The result may not be reliable.',''));
    end
    if sce1.NumGenes < 50 || sce2.NumGenes < 50
        uiwait(warndlg('One of groups contains too few genes (n < 50). The result may not be reliable.',''));
    end

    [T1, ~, ~, xyz1] = sc_splinefit(sce1.X, sce1.g, true, false);
    [T2, ~, ~, xyz2] = sc_splinefit(sce2.X, sce2.g, true, false);
    
    [~,ia,ib] = intersect(T1.genes, T2.genes);
    T1 = T1(ia, :);
    T2 = T2(ib, :);

    %x1 = T1.lgu; y1 = T1.lgcv; z1 = T1.dropr;
    %x2 = T2.lgu; y2 = T2.lgcv; z2 = T2.dropr;
    
    valididx = T1.nearidx>0 & T2.nearidx>0;
    T1 = T1(valididx,:);
    T2 = T2(valididx,:);

    assert(isequal(T1.genes, T2.genes))
    
    bd = vecnorm(xyz1(T1.nearidx,:) - xyz2(T2.nearidx,:),2,2); % baseline difference
    dd = abs(T1.d - T2.d);
    ddn = dd./bd;
    
    BaselineDiffDist = zeros(size(T1,1), 1);
    DiffDistRaw =  zeros(size(T1,1), 1);
    DiffDistNormlized = zeros(size(T1,1), 1);
    
    BaselineDiffDist(valididx) = bd;
    DiffDistRaw(valididx) = dd;
    DiffDistNormlized(valididx) = ddn;
    
    T1.Properties.VariableNames = append(T1.Properties.VariableNames, sprintf('_%s',cL1{1}));
    T2.Properties.VariableNames = append(T2.Properties.VariableNames, sprintf('_%s',cL2{1}));
    
    T=[T1 T2 table(BaselineDiffDist, DiffDistRaw, DiffDistNormlized)];
    
    %assignin("base","T",T)
    assignin("base","T1",T1)
    assignin("base","T2",T2)
    assignin("base","xyz1",xyz1)
    assignin("base","xyz2",xyz2)
    
    T = sortrows(T,"DiffDistNormlized","descend");
    
    %    gui.gui_waitbar(fw);

    outfile = sprintf('%s_vs_%s', ...
        matlab.lang.makeValidName(string(cL1)), matlab.lang.makeValidName(string(cL2)));
    [~, filesaved] = gui.i_exporttable(T, true, 'Tdvgenelist', outfile);
    if ~isempty(filesaved)
       waitfor(helpdlg(sprintf('Result has been saved in %s',filesaved),''));
    end

% ---------

%{
hFig = figure('Visible','off');
hFig.Position(3)=hFig.Position(3)*1.8;
[px_new] = gui.i_getchildpos(FigureHandle, hFig);
if ~isempty(px_new)
    movegui(hFig, px_new);
else
    movegui(hFig, 'center');
end
tb = findall(hFig, 'Tag', 'FigureToolBar'); % get the figure's toolbar handle
uipushtool(tb, 'Separator', 'off');
hFig.Visible=true;
hAx0 = subplot(2,2,[1 3]);
%h1 = scatter3(hAx0, x1, y1, z1, 'filled', 'MarkerFaceAlpha', .1);
hold on
plot3(hAx0, xyz1(:, 1), xyz1(:, 2), xyz1(:, 3), '-', 'linewidth', 4);
%h2 = scatter3(hAx0, x2, y2, z2, 'filled', 'MarkerFaceAlpha', .1);
plot3(hAx0, xyz2(:, 1), xyz2(:, 2), xyz2(:, 3), '-', 'linewidth', 4);

xlabel(hAx0,'Mean+1, log');
ylabel(hAx0,'CV+1, log');
zlabel(hAx0,'Dropout rate (% of zeros)');


    % if ~isempty(g)
    %     dt = datacursormode(hFig);
    %     % dt.UpdateFcn = {@i_myupdatefcn1, g};
    % else
    %     dt = [];
    % end

hAx1 = subplot(2,2,2);
x1 = X1(1,:);
sh1 = plot(hAx1, 1:length(x1), x1);
xlim(hAx1,[1 size(X1,2)]);
title(hAx1, '');
[titxt] = gui.i_getsubtitle(x1);
subtitle(hAx1, titxt);
xlabel(hAx1,'Cell Index');
ylabel(hAx1,'Expression Level');

hAx2 = subplot(2,2,4);
x2 = X2(1,:);
sh2 = plot(hAx2, 1:length(x2), x2);
xlim(hAx2,[1 size(X2,2)]);
title(hAx2, '');
[titxt] = gui.i_getsubtitle(x1);
subtitle(hAx2, titxt);
xlabel(hAx2,'Cell Index');
ylabel(hAx2,'Expression Level');
%}
