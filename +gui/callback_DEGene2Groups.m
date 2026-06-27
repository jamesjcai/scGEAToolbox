function callback_DEGene2Groups(src, ~)

[FigureHandle, sce] = gui.gui_getfigsce(src);
if ~isempty(FigureHandle)
    figure(FigureHandle);
    cleanupObj = onCleanup(@() figure(FigureHandle));
end

extprogname = 'scgeatool_DEAnalysis';
preftagname = 'externalwrkpath';
[wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname, FigureHandle);
if isempty(wkdir), return; end
outdir = wkdir;

[i1, i2, cL1, cL2] = gui.i_select2smplgrps(sce, false, FigureHandle);
if isscalar(i1) || isscalar(i2), return; end

[i1, i2, cL1, cL2, cancelled] = gui.i_whichvswhich(FigureHandle, i1, i2, cL1, cL2);
if cancelled, return; end

fw = gui.myWaitbar(FigureHandle, [], false, 'Computing DE results...');
cleanupFw = onCleanup(@() i_closewaitbar(fw));

try
    T = sc_deg(sce.X(:, i1), sce.X(:, i2), sce.g, 1, true, FigureHandle);
catch ME
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    return;
end

outfile = sprintf('%s_vs_%s_DE_results', ...
    matlab.lang.makeValidName(string(cL1)), ...
    matlab.lang.makeValidName(string(cL2)));
filesaved = fullfile(outdir, [outfile, '.xlsx']);

degparamtag = 'degtestparamset';
paramset = getpref('scgeatoolbox', degparamtag, {0.05, 1.0, 0.01, 'Adjusted P-value'});
[Tup, Tdn, paramset] = pkg.e_processdetable(T, paramset, FigureHandle);
[T, Tnt] = pkg.in_DETableProcess(T, cL1, cL2, sum(i1), sum(i2));

gui.myWaitbar(FigureHandle, fw, false, '', 'Saving DE results...', 0.85);
try
    writetable(T, filesaved, 'FileType', 'spreadsheet', 'Sheet', 'All_processed');
    writetable(Tup, filesaved, 'FileType', 'spreadsheet', 'Sheet', 'Up-regulated');
    writetable(Tdn, filesaved, 'FileType', 'spreadsheet', 'Sheet', 'Down-regulated');
    writetable(Tnt, filesaved, 'FileType', 'spreadsheet', 'Sheet', 'Note');
catch ME
    warning(ME.message);
end

i_closewaitbar(fw);

plotAction(1) = struct('Text', 'Generate Volcano Plot', ...
    'Tooltip', 'Generate volcano plot for DE results', ...
    'Callback', @in_callback_generatevolcano);
plotAction(2) = struct('Text', 'Enrichr Analysis', ...
    'Tooltip', 'Run Enrichr with up/down-regulated DE genes', ...
    'Callback', @in_callback_enrichr_fromtable);
gui.TableViewerApp(T, FigureHandle, outfile, plotAction);


function in_callback_generatevolcano(~, figtab)
    e_volcano(T, Tup, Tdn, figtab);
end

function in_callback_enrichr_fromtable(~, figtab)
    fw2 = gui.myWaitbar(figtab, [], false, 'Running Enrichr analysis...');
    try
        gui.e_enrichrxlsx(Tup, Tdn, T, filesaved);
    catch ME
        gui.myWaitbar(figtab, fw2, true);
        gui.myErrordlg(figtab, ME.message, ME.identifier);
        return;
    end
    gui.myWaitbar(figtab, fw2, true);
    winopen(fileparts(filesaved));
end

function i_closewaitbar(fwx)
    if nargin < 1 || isempty(fwx), return; end
    try
        if pkg.i_isvalid(fwx), close(fwx); end
    catch
        % waitbar may already be closed; safe to ignore
    end
end

function in_callback_savetable(srcx, ~)
    hFig = srcx.Parent.Parent;
    gui.i_exporttable(T, true, ...
        'Tdegenelist', outfile, [], "All_genes", hFig);
end

function in_callback_runenrichr(srcx, ~)
    hFig = srcx.Parent.Parent;
    disp('To run enrichment analysis, type:');
    disp('run.web_Enrichr(Tup.gene(1:250))');
    disp('run.web_Enrichr(Tdn.gene(1:250))');

    [outgenelist, outbackgroundlist, enrichrtype] = ...
        gui.gui_prepenrichr(Tup.gene, sce.g,...
       'Run enrichment analysis with up-regulated DE genes?', ...
       hFig);

    if ~isempty(outbackgroundlist)
        gui.callback_RunEnrichr(src, [], outgenelist, ...
            enrichrtype, ...
            outbackgroundlist, "Up", outdir);
    end

    [outgenelist, outbackgroundlist, enrichrtype] = ...
        gui.gui_prepenrichr(Tdn.gene, sce.g,...
       'Run enrichment analysis with down-regulated DE genes?', ...
       hFig);

    if ~isempty(outbackgroundlist)
        gui.callback_RunEnrichr(src, [], outgenelist, enrichrtype, ...
            outbackgroundlist, "Down", outdir);
    end
end

% ------- begining of volcano_plot

function hFig = e_volcano(T, Tup, Tdn, parentfig)
    T=T(~ismember(T.gene, [Tup.gene; Tdn.gene]),:);
    hx = gui.myFigure(parentfig);
    hFig = hx.FigHandle;
    ax = hx.AxHandle;
    if isempty(ax), ax = gca; end

    hx.addCustomButton('off', @in_callback_runenrichr, 'www.jpg', 'Run enrichment analysis');
    hx.addCustomButton('off', @in_callback_savetable, 'floppy-disk-arrow-in.jpg', 'Export DE Gene Table...');

    hold(ax, "on");
    e_v(Tdn, ax);
    h = e_v(T, ax);
    e_v(Tup, ax);

    xl = max(abs(xlim(ax)));
    xlim(ax, [-xl, xl]);

    h.MarkerFaceColor=[.5 .5 .5];
    h.MarkerEdgeColor=[.5 .5 .5];
    title(ax, sprintf('%s vs. %s', ...
        string(cL1), ...
        string(cL2)));
    ylabel(ax, '-log_{10}(Adj. P-value)')
    xlabel(ax, 'log_{2}(FC)');
    lgd = legend(ax, {sprintf('Down-regulated (%d)', height(Tdn)), ...
        sprintf('Not Sig. (%d)', height(T)), ...
        sprintf('Up-regulated (%d)', height(Tup))},'Location', ...
        'bestoutside');
    try
        mindiffpct = paramset{1};
        minabsolfc = paramset{2};
        apvaluecut = paramset{3};

        Text_below_legend = sprintf('Dropout Diff. > %d%%\nLog2(FC) > %.2f\nAdj. P-Value < %g', ...
            100*mindiffpct, minabsolfc, apvaluecut);
        txt = text(ax, 0, 0, Text_below_legend, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'top', ...
            'FontSize', 10);
        updateTextbox;
        hFig.AutoResizeChildren = 'off';
        hFig.SizeChangedFcn = @updateTextbox;
    catch ME
        disp(ME.message);
    end
    hold(ax, "off");
    hx.show;

    function updateTextbox(~, ~)
        ax.Units = 'normalized';
        axPos = ax.Position;
        lgd.Units = 'normalized';
        lgdPos = lgd.Position;
        legendCenterNorm = [lgdPos(1) + lgdPos(3)/2, ...
                            lgdPos(2)];
        axesNormX = (legendCenterNorm(1) - axPos(1)) / axPos(3);
        axesNormY = (legendCenterNorm(2) - axPos(2)) / axPos(4);
        xLimits = xlim(ax);
        yLimits = ylim(ax);
        xData = xLimits(1) + axesNormX * (xLimits(2) - xLimits(1));
        yData = yLimits(1) + axesNormY * (yLimits(2) - yLimits(1)) - 0.05 * range(yLimits);
        txt.Position = [xData, yData, 0];
    end

    function h = e_v(T, ax)
        genelist = T.gene;
        pvals = T.p_val_adj;
        fc = T.avg_log2FC;
        h = ix_volcanoplot(fc, pvals, genelist, ax);
    end

    function h = ix_volcanoplot(fc, pvals, genelist, ax)
        % Vocano plot
        pvals(pvals < 1e-100) = 1e-100;
        fc(fc<-999) = -10;
        fc(fc>999) = 10;
        x = fc;
        y = -log10(pvals);
        h = scatter(ax, x, y, 8, "filled");

        if isempty(genelist)
            disp('Empty genelist.');
        else
            h.DataTipTemplate.DataTipRows = dataTipTextRow('', genelist);
        end
        if any(abs(xlim(ax))>=10)
            xlim(ax,[-10 10]);
        end
    end
end

% ------- end of volcano_plot

end
