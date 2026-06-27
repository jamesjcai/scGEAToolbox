function callback_RunCoGAPS(src, ~)
% NMF pattern discovery: R Bayesian CoGAPS [PMID:37828301] or MATLAB NMF.
% Falls back to run.ml_cogaps (pure MATLAB) when R is not configured.

[FigureHandle, sce] = gui.gui_getfigsce(src);
if ~gui.gui_showrefinfo('CoGAPS [PMID:37828301]', FigureHandle), return; end

% Choose backend: R (true Bayesian CoGAPS) or MATLAB (fast NMF fallback).
hasR = ispref('scgeatoolbox', 'rexecutablepath') && ...
    ~isempty(getpref('scgeatoolbox', 'rexecutablepath', []));
if hasR
    backend = gui.myQuestdlg(FigureHandle, 'Choose CoGAPS backend:', '', ...
        {'R (Bayesian CoGAPS)', 'MATLAB (fast NMF)'}, 'R (Bayesian CoGAPS)');
    if isempty(backend), return; end
else
    backend = 'MATLAB (fast NMF)';
end
useR = strcmp(backend, 'R (Bayesian CoGAPS)');

wrkdir = [];
if useR
    extprogname = 'R_CoGAPS';
    preftagname = 'externalwrkpath';
    [wrkdir] = gui.i_getwrkdir(preftagname);
    if isempty(wrkdir)
        [wrkdir] = gui.gui_setprgmwkdir(extprogname, preftagname, FigureHandle);
        if isempty(wrkdir), return; end
    end
end

% Select genes (CoGAPS on all genes is slow; HVGs are recommended).
[k, usehvgs] = gui.i_gethvgnum(sce, FigureHandle);
if isempty(k), return; end
if usehvgs
    [~, Xsub, gsub] = sc_hvg(sce.X, sce.g);
    k = min(k, size(Xsub, 1));
    Xsub = Xsub(1:k, :);
    gsub = gsub(1:k);
else
    Xsub = sce.X;
    gsub = sce.g;
end

nPatterns = gui.i_inputnumk(8, 2, 50, 'Number of patterns', FigureHandle);
if isempty(nPatterns), return; end
if useR
    nIterations = gui.i_inputnumk(1000, 100, 50000, ...
        'Number of MCMC iterations (large values are slow)', FigureHandle);
else
    nIterations = gui.i_inputnumk(1000, 100, 50000, ...
        'Max NMF iterations', FigureHandle);
end
if isempty(nIterations), return; end

isdebug = false;
if useR
    answer = gui.myQuestdlg(FigureHandle, 'Keep temporary working files?');
    switch answer
        case 'Yes'
            isdebug = true;
        case 'No'
            isdebug = false;
        otherwise
            return;
    end
end

fw = gui.myWaitbar(FigureHandle);
try
    if useR
        [A, P, Tmarkers] = run.r_cogaps(Xsub, gsub, nPatterns, ...
            nIterations, true, wrkdir, isdebug);
    else
        [A, P, Tmarkers] = run.ml_cogaps(Xsub, gsub, nPatterns, ...
            nIterations, true);
    end
catch ME
    gui.myWaitbar(FigureHandle, fw, true);
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    rethrow(ME);
end
gui.myWaitbar(FigureHandle, fw);

if isempty(P) || size(P, 1) ~= sce.NumCells
    gui.myErrordlg(FigureHandle, 'CoGAPS error: no valid result returned.');
    return;
end

% Store per-cell pattern weights as cell attributes.
for k = 1:size(P, 2)
    sce.setCellAttribute(sprintf('CoGAPS_P%d', k), P(:, k));
end
gui.myGuidata(FigureHandle, sce, src);

% Visualize a chosen pattern on the current embedding.
patternToShow = gui.i_inputnumk(1, 1, size(P, 2), ...
    'Color cells by which pattern?', FigureHandle);
if ~isempty(patternToShow)
    hx = gui.myFigure(FigureHandle);
    gui.i_gscatter3(sce.s, P(:, patternToShow), 1, 1, hx.AxHandle);
    colorbar(hx.AxHandle);
    title(hx.AxHandle, sprintf('CoGAPS Pattern %d', patternToShow));
    hx.show(FigureHandle);
end

% Export factorization matrices and pattern markers.
if ~(ismcc || isdeployed)
    labels = {'Save gene loadings (A matrix) to variable named:', ...
        'Save cell pattern weights (P matrix) to variable named:', ...
        'Save PatternMarker table to variable named:'};
    vars = {'CoGAPS_A', 'CoGAPS_P', 'CoGAPS_PatternMarkers'};
    values = {A, P, Tmarkers};
    export2wsdlg(labels, vars, values);
else
    gui.i_exporttable(Tmarkers, false, 'CoGAPS_PatternMarkers', ...
        [], [], [], FigureHandle);
end

end
