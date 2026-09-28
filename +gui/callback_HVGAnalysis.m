function callback_HVGAnalysis(src, ~, sce)
%CALLBACK_HVGANALYSIS Rank highly variable genes and explore their profiles.
%
%   Behind Analyze > Gene Statistics > Highly Variable Genes (HVGs). Ranks
%   the genes of SCE (or of the figure's dataset) by one of three methods,
%   and shows either the table or the HVG curve plot, where brushed genes
%   can be sent to Enrichr. Was CALLBACK_ENRICHRHVGS, a name that put the
%   optional last step first.

% The window comes from SRC even when SCE is passed in (the brushed-cells
% path): with FigureHandle left empty, every dialog and the progress bar
% here opened unparented.
if nargin < 3 || isempty(sce)
    [FigureHandle, sce] = gui.gui_getfigsce(src);
else
    FigureHandle = gui.gui_getfigsce(src);
end

% Named the same way as GUI.CALLBACK_DVGENE2GROUPS names them, since
% these are the same three rankers. The analytic curve leads: it is what
% SINGLECELLEXPERIMENT.EMBEDCELLS now selects genes with, so the table
% shown here matches the genes an embedding was built from.
optSpline = 'Splinefit Method [PMID:31697351]';
optAnalytic = 'Analytic Curve (closed-form Spline-DV)';
optBrennecke = 'Brennecke et al. (2013) [PMID:24056876]';

answer = gui.myQuestdlg(FigureHandle, 'Which HVG detecting method to use?', '', ...
    {optAnalytic, optSpline, optBrennecke}, optAnalytic);

if ~any(strcmp(answer, {optAnalytic, optSpline, optBrennecke})), return; end

% One of the two, not the table and then the plot over it, which stacked
% three windows. The plot exports its own ranking from its toolbar.
optTable = 'Table';
optPlot = 'Curve Plot';
shown = gui.myQuestdlg(FigureHandle, ['Show the HVG ranking as a table, ', ...
    'or explore it on the curve plot? The plot can export its ranking ', ...
    'as a table from its toolbar.'], '', ...
    {optTable, optPlot}, optPlot);
if ~any(strcmp(shown, {optTable, optPlot})), return; end

if strcmp(shown, optPlot)
    try
        switch answer
            case optAnalytic
                gui.i_hvgcurveplot(sce.X, sce.g, true, true, ...
                    FigureHandle, "analytic");
            case optSpline
                gui.i_hvgcurveplot(sce.X, sce.g, true, true, ...
                    FigureHandle, "splinefit");
            otherwise % optBrennecke; the other answers returned above
                sc_hvg(sce.X, sce.g, true, true);
        end
    catch ME
        gui.myErrordlg(FigureHandle, ME.message);
    end
    return;
end

fw = gui.myWaitbar(FigureHandle);
try
    switch answer
        case optAnalytic
            T = sc_analyticfit(sce.X, sce.g);
        case optSpline
            T = sc_splinefit(sce.X, sce.g, true, false);
        otherwise % optBrennecke; the other answers returned above
            T = sc_hvg(sce.X, sce.g, true, false);
    end
catch ME
    gui.myWaitbar(FigureHandle, fw, true);
    gui.myErrordlg(FigureHandle, ME.message);
    return;
end
gui.myWaitbar(FigureHandle, fw);
gui.TableViewerApp(T, FigureHandle);

end
