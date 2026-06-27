function callback_DVGene2Groups(src, ~)

lcolors = lines(2);
lcolor1 = lcolors(1,:);
lcolor2 = lcolors(2,:);

[FigureHandle, sce_ori] = gui.gui_getfigsce(src);
sce = copy(sce_ori);

if ~gui.gui_showrefinfo('DV Analysis', FigureHandle), return; end

extprogname = 'scgeatool_DVAnalysis';
preftagname = 'externalwrkpath';
[wrkdir] = gui.gui_setprgmwkdir(extprogname, preftagname, FigureHandle);
if isempty(wrkdir), return; end


a=sce.NumGenes;
[sce] = gui.i_selectinfogenes(sce, [], FigureHandle);
b=sce.NumGenes;
fprintf('%d genes removed.\n', a-b);

[i1, i2, cL1, cL2] = gui.i_select2smplgrps(sce, false, FigureHandle);
if isscalar(i1) || isscalar(i2), return; end


c=zeros(size(i1));
c(i1)=1; c(i2)=2;

cL1 = matlab.lang.makeValidName(cL1);
cL2 = matlab.lang.makeValidName(cL2);
    if isequal(cL1, cL2)
        tmptxtc = matlab.lang.makeUniqueStrings([cL1 cL2]);
        cL1 = tmptxtc(1);
        cL2 = tmptxtc(2);
    end
    if ~all(c>0)
        sce = sce.selectcells(c>0); % OK
        c=c(c>0);
    end
i1=c==1;
i2=c==2;

sce1 = copy(sce);
sce1 = sce1.selectcells(i1); % OK
sce1 = sce1.qcfilter; % OK


sce2 = copy(sce);
sce2 = sce2.selectcells(i2); % OK
sce2 = sce2.qcfilter; % OK


    if sce1.NumCells < 10 || sce2.NumCells < 10 || sce1.NumGenes < 10 || sce2.NumGenes < 10
        gui.myErrordlg(FigureHandle, ['Filtered SCE contains too' ...
            ' few cells (n < 10) or genes (n < 10).'],'','modal');
        return;
    end

    % assignin('base', "sce1", sce1);
    % assignin('base', "sce2", sce2);
    % assignin('base', "cL1", cL1);
    % assignin('base', "cL2", cL2);

a = 'Splinefit Method [PMID:40113778]';
b = 'Brennecke et al. (2013) [PMID:24056876]';

            % answerx = gui.myQuestdlg(FigureHandle, ...
            %     'Which HVG detecting method to use?', '', ...
            %     {a, b}, a);
answerx = a;

fw = gui.myWaitbar(FigureHandle, [], false, 'Computing DV results...');
cleanupObj = onCleanup(@() i_closewaitbar(fw)); 

try
    switch answerx
        case a
            [T, X1, X2, g, xyz1, xyz2, ...
                px1, py1, pz1, ...
                px2, py2, pz2] = sc_dvg(sce1, sce2, cL1, cL2, 'splinefit');
            methodtag = 'splinefit';
        case b
            T = sc_dvg(sce1, sce2, cL1, cL2, 'brennecke');
            methodtag = 'brennecke';
        otherwise
            return;
    end
catch ME
    gui.myWaitbar(FigureHandle, fw, true);
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    return;
end

outfile = sprintf('%s_vs_%s_DV_%s_results', ...
        matlab.lang.makeValidName(string(cL1)), ...
        matlab.lang.makeValidName(string(cL2)), ...
        methodtag);
filesaved = fullfile(wrkdir, [outfile, '.xlsx']);

drawnow;
gui.myWaitbar(FigureHandle, fw, false, '', 'Saving DV results...', 0.85);
try
    writetable(T, filesaved, 'FileType', 'spreadsheet', 'Sheet', 'DV_results');
catch ME
    warning(ME.message);
end
i_closewaitbar(fw);

hFig = []; hx = []; h1 = []; h2 = [];
hAx0 = []; hAx1 = []; hAx2 = [];
sh1 = []; sh2 = []; idx = 1;
h3 = []; h3a = []; h3b = []; h4 = []; h5 = [];

enrichrAction = struct('Text', 'Enrichr Analysis', ...
    'Tooltip', 'Run Enrichr with DV genes', ...
    'Callback', @in_callback_enrichr_fromtable);
if strcmp(answerx, a)
    plotAction(1) = struct('Text', 'Open Plot', ...
        'Tooltip', 'Open the interactive DV scatter plot', ...
        'Callback', @in_callback_openDVplot);
    plotAction(2) = enrichrAction;
    figtab = gui.TableViewerApp(T, FigureHandle, outfile, plotAction);
else
    figtab = gui.TableViewerApp(T, FigureHandle, outfile, enrichrAction);
end


function in_callback_openDVplot(~, ~)
        hx = gui.myFigure(figtab, true);
        hFig = hx.FigHandle;
        hFig.Position(3) = hFig.Position(3)*1.8;

        hx.addCustomButton('off', {@in_callback_HighlightSelectedGenes, 1}, 'list.gif', 'Selet a gene to show expression profile');
        hx.addCustomButton('off', {@in_callback_HighlightSelectedGenes, 2}, 'list2.gif', 'Selet a gene from sorted list');
        hx.addCustomButton('off', @in_callback_viewTable, 'icon-fa-stack-exchange-10.gif', 'View DV gene table...');
        hx.addCustomButton('on', @in_callback_EnrichrHVGs, 'plotpicker-andrewsplot.gif', 'Select top n genes to perform web-based enrichment analysis...');
        hx.addCustomButton('off', @in_callback_Enrichr, 'plotpicker-andrewsplot.gif', 'Enrichr test...');
        hx.addCustomButton('off', @in_callback_genecards, 'www.jpg', 'GeneCards...');
        hx.addCustomButton('on', @in_callback_ChangeAlphaValue, 'plotpicker-rose.gif', 'Change MarkerFaceAlpha value');
        hx.addCustomButton('off', @in_callback_changeMarkerSize, 'icon-mat-text-fields-10.gif', 'ChangeFontSize');

        hAx0 = subplot(2,2,[1 3]);
        h1 = scatter3(hAx0, px1, py1, pz1, 'filled', 'MarkerFaceAlpha', .1);
        hold on
        h2 = scatter3(hAx0, px2, py2, pz2, 'filled', 'MarkerFaceAlpha', .1);
        plot3(hAx0, xyz1(:, 1), xyz1(:, 2), xyz1(:, 3), '-', 'linewidth', 4, 'Color', lcolor1);
        plot3(hAx0, xyz2(:, 1), xyz2(:, 2), xyz2(:, 3), '-', 'linewidth', 4, 'Color', lcolor2);

        xlabel(hAx0, 'Mean, log1p');
        ylabel(hAx0, 'CV, log1p');
        zlabel(hAx0, 'Dropout rate (% of zeros)');

        if ~isempty(g)
            dt = datacursormode(hFig);
            datacursormode(hFig, 'on');
            dt.UpdateFcn = {@in_myupdatefcn3, g};
        end

        idx = find(g == table2array(T(1, 1)));

        hAx1 = subplot(2,2,2);
        x1 = X1(idx,:);
        sh1 = plot(hAx1, 1:length(x1), x1, 'Color', lcolor1);
        xlim(hAx1, [1 size(X1,2)]);

        title(hAx1, strrep(sprintf('%s', g(idx)), '_', '\_'));
        [titxt] = gui.i_getsubtitle(x1, cL1{1});
        subtitle(hAx1, titxt);
        xlabel(hAx1, 'Cell Index');
        ylabel(hAx1, 'Expression Level');

        hAx2 = subplot(2,2,4);
        x2 = X2(idx,:);
        sh2 = plot(hAx2, 1:length(x2), x2, 'Color', lcolor2);
        xlim(hAx2, [1 size(X2,2)]);

        title(hAx2, strrep(sprintf('%s', g(idx)), '_', '\_'));
        [titxt] = gui.i_getsubtitle(x2, cL2{1});
        subtitle(hAx2, titxt);
        xlabel(hAx2, 'Cell Index');
        ylabel(hAx2, 'Expression Level');
        h3 = []; h4 = []; h5 = [];
        h3a = []; h3b = [];
        yl = cell2mat(get([hAx1, hAx2], 'Ylim'));
        ylnew = [min(yl(:, 1)), max(yl(:, 2))];
        set([hAx1, hAx2], 'Ylim', ylnew);
        subplot(hAx0);
        hx.show(figtab);
    end

function in_callback_Enrichr(~, ~)
        answer = gui.myQuestdlg(hFig, 'Enrichr test with top DV genes. Continue?','');
        if ~strcmp(answer,'Yes'), return; end
        answer = gui.myQuestdlg(hFig, 'Select type of DV genes.','',...
            {'Mixed','Varibility increasing','Varibility decreasing'},'Mixed');
        switch answer
            case 'Mixed'
               Tin = T;
            case 'Varibility increasing'
                Tup = T(T.DiffSign > 0, :);
                Tin = Tup;
            case 'Varibility decreasing'
                Tdn = T(T.DiffSign < 0, :);
                Tin = Tdn;
        end

        [outgenelist, outbackgroundlist, enrichrtype] = ...
            gui.gui_prepenrichr(Tin.gene(1:250), Tin.gene,...
                sprintf('Run enrichment analysis with %s DV genes?', lower(answer)), ...
                hFig);
        gui.callback_RunEnrichr(src, [], outgenelist, enrichrtype, outbackgroundlist);

    end

function in_callback_viewTable(~, ~)
        if ~isempty(figtab) && pkg.i_isvalid(figtab)
            gui.i_bringtofront(figtab);
        else
            figtab = gui.TableViewerApp(T, hx.FigHandle, outfile);
        end
    end

function txt = in_myupdatefcn3(src, event_obj, g)
        if isequal(get(src, 'Parent'), hAx0)
            subplot(hAx0);
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

            x_cleanfigspace(false);

            % if ~isempty(h3), delete(h3); end
            % if ~isempty(h3a), delete(h3a); end
            % if ~isempty(h3b), delete(h3b); end
            % if ~isempty(h4), delete(h4); end
            % if ~isempty(h5), delete(h5); end

                h3 = plot3(hAx0, [px1(idx) px2(idx)], ...
                    [py1(idx), py2(idx)], ...
                    [pz1(idx), pz2(idx)],'k-','LineWidth',1);

                h3a = plot3(hAx0, px1(idx), py1(idx), pz1(idx), 'b.','MarkerSize',10);
                h3b = plot3(hAx0, px2(idx), py2(idx), pz2(idx), 'r.','MarkerSize',10);
                % daspect(hAx0,'auto');
                % h3b = arrow3([px1(idx), py1(idx), pz1(idx)], ...
                %    [px2(idx), py2(idx), pz2(idx)]);
                % h3b = quiver3(hAx0, px1(idx), py1(idx), pz1(idx), px2(idx), py2(idx), pz2(idx));
            txt = {g(idx)};
            x1 = X1(idx, :);
            if ~isempty(sh1) && pkg.i_isvalid(sh1), delete(sh1); end
            sh1 = plot(hAx1, 1:length(x1), x1, 'Color',lcolor1);
            xlim(hAx1,[1 size(X1,2)]);
            title(hAx1, strrep(sprintf('%s',g(idx)),'_','\_') );
            [titxt] = gui.i_getsubtitle(x1, cL1{1});
            subtitle(hAx1, titxt);
            xlabel(hAx1,'Cell Index');
            ylabel(hAx1,'Expression Level');

            x2 = X2(idx, :);
            if ~isempty(sh2) && pkg.i_isvalid(sh2), delete(sh2); end
            sh1 = plot(hAx2, 1:length(x2), x2, 'Color', lcolor2);
            xlim(hAx2,[1 size(X2,2)]);
            title(hAx2, strrep(sprintf('%s',g(idx)),'_','\_'));
            [titxt] = gui.i_getsubtitle(x2, cL2{1});
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

function x_cleanfigspace(deldatatip)
        if nargin<1, deldatatip = false; end
        if deldatatip
           dtp = findobj(h1, 'Type', 'datatip');
           if ~isempty(dtp), delete(dtp); end
           dtp = findobj(h2, 'Type', 'datatip');
           if ~isempty(dtp), delete(dtp); end
        end

       if ~isempty(h3), delete(h3); end
       if ~isempty(h3a), delete(h3a); end
       if ~isempty(h3b), delete(h3b); end
       if ~isempty(h4), delete(h4); end
       if ~isempty(h5), delete(h5); end
    end

function in_callback_HighlightSelectedGenes(~, ~, typeid)
       if nargin < 3, typeid = 1; end

       x_cleanfigspace(true);

       % dtp = findobj(h1, 'Type', 'datatip');
       % if ~isempty(dtp), delete(dtp); end
       % dtp = findobj(h2, 'Type', 'datatip');
       % if ~isempty(dtp), delete(dtp); end
       % if ~isempty(h3), delete(h3); end
       % if ~isempty(h3a), delete(h3a); end
       % if ~isempty(h3b), delete(h3b); end
       % if ~isempty(h4), delete(h4); end
       % if ~isempty(h5), delete(h5); end

       switch typeid
           case 1
                gsorted = natsort(g);
           case 2
                gsorted = T.(T.Properties.VariableNames{1});
       end
        if gui.i_isuifig(FigureHandle)
            [indx2, tf2] = gui.myListdlg(FigureHandle, gsorted, ...
                'Select a gene:');
        else
            [indx2, tf2] = listdlg('PromptString', ...
                'Select a gene:', ...
                'SelectionMode', 'single', ...
                'ListString', gsorted, ...
                'ListSize', [220, 300]);
        end
        if tf2 == 1
            idx = find(g == gsorted{indx2});
        else
            return;
        end
        h1.BrushData = idx;
        h2.BrushData = idx;
        datatip(h1, 'DataIndex', idx);
        datatip(h2, 'DataIndex', idx);

        x_cleanfigspace(false);
        % if ~isempty(h3), delete(h3); end
        % if ~isempty(h4), delete(h4); end
        % if ~isempty(h5), delete(h5); end

         % h3 = plot3(hAx0, [px1(idx) px2(idx)], ...
         %     [py1(idx), py2(idx)], ...
         %     [pz1(idx), pz2(idx)],'k--','LineWidth',1);

         [nearidx] = dsearchn(xyz1, [px1(idx) py1(idx) pz1(idx)]);

         % assignin('base',"xyz1",xyz1);
         % h4 = arrow3(xyz1(nearidx, :), [px1(idx), py1(idx), pz1(idx)]);

        h4 = plot3(hAx0, [px1(idx) xyz1(nearidx, 1)], ...
            [py1(idx), xyz1(nearidx, 2)], ...
            [pz1(idx), xyz1(nearidx, 3)],'-','LineWidth',2,'Color',lcolor1);

         [nearidx] = dsearchn(xyz2, [px2(idx) py2(idx) pz2(idx)]);

         % h5 = arrow3(xyz2(nearidx, :), [px2(idx), py2(idx), pz2(idx)]);

         h5 = plot3(hAx0, [px2(idx) xyz2(nearidx, 1)], ...
             [py2(idx), xyz2(nearidx, 2)], ...
             [pz2(idx), xyz2(nearidx, 3)],'k-','LineWidth',2,'Color',lcolor2);

    end

function in_callback_EnrichrHVGs(~, ~)
        k = gui.i_inputnumk(200, 1, 2000, 'Select top n genes', FigureHandle);
        if ~isempty(k)
            gsorted = T.(T.Properties.VariableNames{1});
            gselected = gsorted(1:k);
            fprintf('%d genes are selected.\n', length(gselected));
            gui.i_enrichtest(gselected, gsorted, k);
        end
    end

function in_callback_genecards(~, ~)
        web(sprintf('https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s', g(idx)),'-new');
    end

function in_callback_changeMarkerSize(~, ~)
        % h1.Marker
        if h1.SizeData > 40
            h1.SizeData = 10;
            h2.SizeData = 10;
        else
            h1.SizeData = h1.SizeData + 2;
            h2.SizeData = h2.SizeData + 2;
        end
    end

function in_callback_ChangeAlphaValue(~, ~)
        if h1.MarkerFaceAlpha <= 0.05
            h1.MarkerFaceAlpha = 1;
            h2.MarkerFaceAlpha = 1;
        else
            h1.MarkerFaceAlpha = h1.MarkerFaceAlpha - 0.1;
            h2.MarkerFaceAlpha = h2.MarkerFaceAlpha - 0.1;
        end
    end

function in_callback_enrichr_fromtable(~, figtab)
    libs = ["GO_Biological_Process_2025", "GO_Molecular_Function_2025", ...
            "KEGG_2021_Human", "Reactome_Pathways_2024"];
    suffixes = ["GO_BP", "GO_MF", "KEGG", "Reactome"];
    groups = { ...
        T,                      "Mix"; ...
        T(T.DiffSign > 0, :),   "VInc"; ...
        T(T.DiffSign < 0, :),   "VDec"};

    fw2 = gui.myWaitbar(figtab, [], false, 'Running Enrichr analysis...');
    try
        for kg = 1:size(groups, 1)
            Tg = groups{kg, 1};
            prefix = groups{kg, 2};
            if height(Tg) == 0, continue; end
            n = min(250, height(Tg));
            Tlist = run.ml_Enrichr(Tg.gene(1:n), T.gene, libs);
            for kl = 1:numel(Tlist)
                if ~isempty(Tlist{kl}) && istable(Tlist{kl}) && height(Tlist{kl}) > 0
                    writetable(Tlist{kl}, filesaved, ...
                        'FileType', 'spreadsheet', ...
                        'Sheet', sprintf('%s_%s', prefix, suffixes(kl)));
                end
            end
        end
    catch ME
        gui.myWaitbar(figtab, fw2, true);
        gui.myErrordlg(figtab, ME.message, ME.identifier);
        return;
    end
    gui.myWaitbar(figtab, fw2, true);
    winopen(fileparts(filesaved));
end

end

function i_closewaitbar(fw)
if nargin < 1 || isempty(fw)
    return;
end

try
    if pkg.i_isvalid(fw)
        close(fw);
    end
catch
    % waitbar may already be closed; safe to ignore
end
end
