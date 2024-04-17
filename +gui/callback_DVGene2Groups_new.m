function callback_DVGene2Groups_new(src, ~)

cx=lines(2);
cx1=cx(1,:);
cx2=cx(2,:);

    FigureHandle = src.Parent.Parent;
    sce = guidata(FigureHandle);
    if ~gui.gui_showrefinfo('DV Analysis'), return; end

% ---------------------------------
spciestag = gui.i_selectspecies(2);
if isempty(spciestag), return; end

prompt = {'Remove Mt-Genes (MT-ND1, MT-ND6, MT-CYB, MT-COI, MT-ATP6, etc.)?', ...
    'Remove Hemoglobin Genes (HBA1, HBB, Hba-a1, etc.)?', ...
    'Remove Ribosomal Genes (RPSA, RPS2, RPS3, RPL3, RPL4, RPLP1, etc.)?', ...
    'Remove Genes Without Approved Symbols?'};
dlgtitle = '';
dims = [1, 85];

definput = {'Yes', 'Yes', 'Yes', 'Yes'};
answer = inputdlg(prompt, dlgtitle, dims, definput);
if isempty(answer), return; end

if strcmpi(answer{1},'Yes') || strcmpi(answer{1},'Y')
    sce = sce.rmmtgenes;
    disp('Mt-genes removed.');    
end

if strcmpi(answer{2},'Yes') || strcmpi(answer{2},'Y')
    sce = sce.rmhemoglobingenes;
    disp('Hemoglobin genes removed.');    
end

if strcmpi(answer{3},'Yes') || strcmpi(answer{3},'Y')
    sce = sce.rmribosomalgenes;
    disp('Ribosomal genes removed.');    
end

if strcmpi(answer{4},'Yes') || strcmpi(answer{4},'Y')
    mfolder = fileparts(mfilename('fullpath'));
    switch spciestag
        case 'human'
            load(fullfile(mfolder, ...
                '../resources', 'Biomart_human_genes.mat'), 'T');
        case 'mouse'
            load(fullfile(mfolder, ...
                '../resources', 'Biomart_mouse_genes.mat'), 'T');
    end
    ApprovedSymbol = string(T.GeneName);
    [idx] = ismember(upper(sce.g), upper(ApprovedSymbol));
    a1 = length(sce.g);
    sce.g(~idx) = [];
    sce.X(~idx, :) = [];
    a2 = length(sce.g);
    fprintf('%d genes without approved symbols are found and removed.\n',a1-a2);
end

% ---------------------------------
    
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

    if ~isequal(sce1.g, sce2.g)
        [g, ia, ib] = intersect(sce1.g, sce2.g,'stable');
        X1 = sce1.X(ia, :);
        X2 = sce2.X(ib, :);
    else
        g = sce1.g;
        X1 = sce1.X;
        X2 = sce2.X;        
    end
    
    X1 = sc_norm(X1,'type','libsize');
    X2 = sc_norm(X2,'type','libsize');

    [T1, X1, g1, xyz1] = sc_splinefit(X1, g, true, false);
    [T1, idx] = sortrows(T1,'genes','descend');
    X1 = X1(idx, :);
    g1 = g1(idx);

    [T2, X2, g2, xyz2] = sc_splinefit(X2, g, true, false);
    [T2, idx] = sortrows(T2,'genes','descend');
    X2 = X2(idx, :);
    g2 = g2(idx);

    assert(isequal(g1, g2))
    g = g1;


    px1 = T1.lgu; py1 = T1.lgcv; pz1 = T1.dropr;
    px2 = T2.lgu; py2 = T2.lgcv; pz2 = T2.dropr;

    d1=([px1 py1 pz1] - xyz1(T1.nearidx,:));
    d2=([px2 py2 pz2] - xyz2(T2.nearidx,:));

    DiffDist = vecnorm(d1 - d2, 2, 2);    
   
    T1.Properties.VariableNames = append(T1.Properties.VariableNames, sprintf('_%s',cL1{1}));
    T2.Properties.VariableNames = append(T2.Properties.VariableNames, sprintf('_%s',cL2{1}));
    
    T = [T1 T2 table(DiffDist)];    
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
% tb = findall(hFig, 'Tag', 'FigureToolBar'); % get the figure's toolbar handle
% uipushtool(tb, 'Separator', 'off');

pkg.i_addbutton2fig(tb, 'off', {@in_HighlightGenes, 1}, 'list.gif', 'Selet a gene to show expression profile');
% pkg.i_addbutton2fig(tb, 'off', @in_HighlightSelectedGenes, 'xplotpicker-qqplot.gif', 'Highlight selected genes');
pkg.i_addbutton2fig(tb, 'off', {@in_HighlightGenes, 2}, 'plotpicker-qqplot.gif', 'Selet a gene from sorted list');
pkg.i_addbutton2fig(tb, 'off', @EnrichrHVGs, 'plotpicker-andrewsplot.gif', 'Select top n genes to perform web-based enrichment analysis...');


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
h3a = [];
h3b = [];
yl = cell2mat(get([hAx1, hAx2], 'Ylim'));
ylnew = [min(yl(:, 1)), max(yl(:, 2))];
set([hAx1, hAx2], 'Ylim', ylnew);


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
            if ~isempty(h3a), delete(h3a); end
            if ~isempty(h3b), delete(h3b); end
                h3 = plot3(hAx0, [px1(idx) px2(idx)], ...
                    [py1(idx), py2(idx)], ...
                    [pz1(idx), pz2(idx)],'k-','LineWidth',1);
                h3a = plot3(hAx0, px1(idx), py1(idx), pz1(idx), 'b.','MarkerSize',10);
                h3b = plot3(hAx0, px2(idx), py2(idx), pz2(idx), 'r.','MarkerSize',10);
                % daspect(hAx0,'auto');
                %h3b = arrow3([px1(idx), py1(idx), pz1(idx)], ...
                %    [px2(idx), py2(idx), pz2(idx)]);
                %h3b = quiver3(hAx0, px1(idx), py1(idx), pz1(idx), px2(idx), py2(idx), pz2(idx));
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

            yl = cell2mat(get([hAx1, hAx2], 'Ylim'));
            ylnew = [min(yl(:, 1)), max(yl(:, 2))];
            set([hAx1, hAx2], 'Ylim', ylnew);
            
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
       if ~isempty(h3a), delete(h3a); end
       if ~isempty(h3b), delete(h3b); end
       
       switch typeid
           case 1
                gsorted = natsort(g);
           case 2
                gsorted = T.(T.Properties.VariableNames{1});
       end

        [indx2, tf2] = listdlg('PromptString', ...
            'Select a gene:', ...
            'SelectionMode', 'single', ...
            'ListString', gsorted, ...
            'ListSize', [220, 300]);
        if tf2 == 1
            idx = find(g == gsorted{indx2});
        else
            return;
        end        
        h1.BrushData = idx;
        h2.BrushData = idx;
        datatip(h1, 'DataIndex', idx);
        datatip(h2, 'DataIndex', idx);
         if ~isempty(h3), delete(h3); end
         h3 = plot3(hAx0, [px1(idx) px2(idx)], ...
             [py1(idx), py2(idx)], ...
             [pz1(idx), pz2(idx)],'k-','LineWidth',1);       
    end

    function EnrichrHVGs(~, ~)
        k = gui.i_inputnumk(150, 1, 2000, 'Select top n genes');
        if ~isempty(k)
            gsorted = T.(T.Properties.VariableNames{1});
            gselected = gsorted(1:k);
            fprintf('%d genes are selected.\n', length(gselected));        
            gui.i_enrichtest(gselected, gsorted, k);
        end
    end

end
