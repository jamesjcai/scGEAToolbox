function callback_DVGene2Groups(src, ~)

%warndlg('The function is under development.');
%return;

cx=lines(2);
cx1=cx(1,:);
cx2=cx(2,:);

    FigureHandle = src.Parent.Parent;
    sce = guidata(FigureHandle);
    if ~gui.gui_showrefinfo('DV Analysis'), return; end
    
    [i1, i2, cL1, cL2] = gui.i_select2grps(sce, false);
    if length(i1) == 1 || length(i2) == 1, return; end

    fw = gui.gui_waitbar;
    
    c=zeros(size(i1));
    c(i1)=1; c(i2)=2;
    % cL=[cL1;cL2];
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

    [T1, X1, g1, xyz1] = sc_splinefit(sce1.X, sce1.g, true, false);
    [T2, X2, g2, xyz2] = sc_splinefit(sce2.X, sce2.g, true, false);
    
    [g, ia, ib] = intersect(T1.genes, T2.genes);
    T1 = T1(ia, :);
    T2 = T2(ib, :);

    X1 = X1(ia, :);
    X2 = X2(ib, :);
    g1 = g1(ia);
    g2 = g2(ib);
 
    valididx = T1.nearidx>0 & T2.nearidx>0;
    T1 = T1(valididx,:);
    T2 = T2(valididx,:);

    X1 = X1(valididx, :);
    g1 = g1(valididx);
    X2 = X2(valididx, :);
    g2 = g2(valididx);
    g = g(valididx);

    assert(isequal(g1, g2))
    assert(isequal(g, g2))
    assert(isequal(T1.genes, T2.genes))
    assert(isequal(T1.genes, g))
    
    %bd = vecnorm(xyz1(T1.nearidx,:) - xyz2(T2.nearidx,:),2,2); % baseline difference
    %dd = abs(T1.d - T2.d);
    %ddn = dd./(1+bd);

    px1 = T1.lgu; py1 = T1.lgcv; pz1 = T1.dropr;
    px2 = T2.lgu; py2 = T2.lgcv; pz2 = T2.dropr;

    d1=([px1 py1 pz1] - xyz1(T1.nearidx,:));
    d2=([px2 py2 pz2] - xyz2(T2.nearidx,:));
    ddn = vecnorm(d1 - d2, 2, 2);


    DiffDist = zeros(size(T1,1), 1);
   
    DiffDist(valididx) = ddn;

    
    T1.Properties.VariableNames = append(T1.Properties.VariableNames, sprintf('_%s',cL1{1}));
    T2.Properties.VariableNames = append(T2.Properties.VariableNames, sprintf('_%s',cL2{1}));
    
    T=[T1 T2 table(DiffDist)];
    
    %assignin("base","T",T)
    %assignin("base","T1",T1)
    %assignin("base","T2",T2)
    %assignin("base","xyz1",xyz1)
    %assignin("base","xyz2",xyz2)
    
    T = sortrows(T,"DiffDist","descend");
    
    gui.gui_waitbar(fw);

    outfile = sprintf('%s_vs_%s', ...
        matlab.lang.makeValidName(string(cL1)), matlab.lang.makeValidName(string(cL2)));
    [~, filesaved] = gui.i_exporttable(T, true, 'Tdvgenelist', outfile);
    if ~isempty(filesaved)
       waitfor(helpdlg(sprintf('Result has been saved in %s',filesaved),''));
    end

% ---------
answer = questdlg('Show gene expression profiles?','');
if ~strcmp(answer,'Yes'), return; end


hFig = figure('Visible','off');
hFig.Position(3) = hFig.Position(3)*1.8;
gui.i_movegui2parent(hFig, FigureHandle);

delete(findall(hFig, 'Tag', 'FigureToolBar'))
tb = uitoolbar('Parent', hFig);
%tb = findall(hFig, 'Tag', 'FigureToolBar'); % get the figure's toolbar handle
uipushtool(tb, 'Separator', 'off');

pkg.i_addbutton2fig(tb, 'off', {@in_HighlightGenes, 1}, 'list.gif', 'Highlight top HVGs');
% pkg.i_addbutton2fig(tb, 'off', @in_HighlightSelectedGenes, 'xplotpicker-qqplot.gif', 'Highlight selected genes');
pkg.i_addbutton2fig(tb, 'off', {@in_HighlightGenes, 2}, 'plotpicker-qqplot.gif', 'Highlight top HVGs');
pkg.i_addbutton2fig(tb, 'off', @EnrichrHVGs, 'plotpicker-andrewsplot.gif', 'Enrichment analysis...');


hFig.Visible=true;
hAx0 = subplot(2,2,[1 3]);
h1 = scatter3(hAx0, px1, py1, pz1, 'filled', 'MarkerFaceAlpha', .1);
hold on
h2 = scatter3(hAx0, px2, py2, pz2, 'filled', 'MarkerFaceAlpha', .1);
plot3(hAx0, xyz1(:, 1), xyz1(:, 2), xyz1(:, 3), '-', 'linewidth', 4, 'Color',cx1);
plot3(hAx0, xyz2(:, 1), xyz2(:, 2), xyz2(:, 3), '-', 'linewidth', 4, 'Color',cx2);

xlabel(hAx0,'Mean+1, log');
ylabel(hAx0,'CV+1, log');
zlabel(hAx0,'Dropout rate (% of zeros)');


if ~isempty(g)
    dt = datacursormode(hFig);
    datacursormode(hFig, 'on');
    dt.UpdateFcn = {@in_myupdatefcn3, g};
end

hAx1 = subplot(2,2,2);
x1 = X1(1,:);
sh1 = plot(hAx1, 1:length(x1), x1, 'Color',cx1);
xlim(hAx1,[1 size(X1,2)]);
title(hAx1, sce1.g(1));
[titxt] = gui.i_getsubtitle(x1);
subtitle(hAx1, titxt);
xlabel(hAx1,'Cell Index');
ylabel(hAx1,'Expression Level');

hAx2 = subplot(2,2,4);
x2 = X2(1,:);
sh2 = plot(hAx2, 1:length(x2), x2, 'Color',cx2);
xlim(hAx2,[1 size(X2,2)]);
title(hAx2, sce2.g(1));
[titxt] = gui.i_getsubtitle(x1);
subtitle(hAx2, titxt);
xlabel(hAx2,'Cell Index');
ylabel(hAx2,'Expression Level');
h3 = [];

    function txt = in_myupdatefcn3(src, event_obj, g)
        if isequal(get(src, 'Parent'), hAx0)

            % dtp = findobj(h1, 'Type', 'datatip');
            % if ~isempty(dtp), delete(dtp); end
            % dtp = findobj(h2, 'Type', 'datatip');
            % if ~isempty(dtp), delete(dtp); end            
            idx = event_obj.DataIndex;
            if idx > length(g)*2
                txt = num2str(event_obj.Position(2)); 
                return; 
            end
            if idx > length(g)
                idx = idx - length(g);
            end
            
            if ~isempty(h3), delete(h3); end
            h3 = plot3(hAx0, [px1(idx) px2(idx)], ...
                [py1(idx), py2(idx)], ...
                [pz1(idx), pz2(idx)],'k-','LineWidth',1);

            txt = {g(idx)};

            x1 = X1(idx, :);
            if ~isempty(sh1) && isvalid(sh1), delete(sh1); end
            sh1 = plot(hAx1, 1:length(x1), x1, 'Color',cx1);
            xlim(hAx1,[1 size(X1,2)]);
            title(hAx1, g(idx));    
            [titxt] = gui.i_getsubtitle(x1);
            subtitle(hAx1, titxt);
            xlabel(hAx1,'Cell Index');
            ylabel(hAx1,'Expression Level');

            x2 = X2(idx, :);
            if ~isempty(sh2) && isvalid(sh2), delete(sh2); end
            sh1 = plot(hAx2, 1:length(x2), x2, 'Color',cx2);
            xlim(hAx2,[1 size(X2,2)]);
            title(hAx2, g(idx));    
            [titxt] = gui.i_getsubtitle(x2);
            subtitle(hAx2, titxt);
            xlabel(hAx2,'Cell Index');
            ylabel(hAx2,'Expression Level');            
        else
            txt = num2str(event_obj.Position(2));
        end
    end

    function in_HighlightGenes(~, ~, typeid)
        if nargin < 3, typeid = 1; end
        %h.MarkerIndices=idx20;
       dtp = findobj(h1, 'Type', 'datatip');
       if ~isempty(dtp), delete(dtp); end
       dtp = findobj(h2, 'Type', 'datatip');
       if ~isempty(dtp), delete(dtp); end
       if ~isempty(h3), delete(h3); end
        
       switch typeid
           case 1
                gsorted = natsort(g);
           case 2
                gsorted = T.(T.Properties.VariableNames{1});
       end

        [indx2, tf2] = listdlg('PromptString', ...
            'Select a gene:', ...
            'SelectionMode', 'single', ...
            'ListString', gsorted);
        if tf2 == 1
            idx = find(g == gsorted{indx2});
        else
            return;
        end        
        h1.BrushData = idx;
        h2.BrushData = idx;
        datatip(h1, 'DataIndex', idx);
        datatip(h2, 'DataIndex', idx);
        % if ~isempty(h3), delete(h3); end
        % h3 = plot3(hAx0, [px1(idx) px2(idx)], ...
        %     [py1(idx), py2(idx)], ...
        %     [pz1(idx), pz2(idx)],'k-','LineWidth',1);       
    end

    function EnrichrHVGs(~, ~)
        k = gui.i_inputnumk(10, 1, 2000);
        if ~isempty(k)
            gsorted = T.(T.Properties.VariableNames{1});
            gselected = gsorted(1:k);
            fprintf('%d genes are selected.\n', length(gselected));        
            gui.i_enrichtest(gselected, gsorted, k);
        end
    end

end
