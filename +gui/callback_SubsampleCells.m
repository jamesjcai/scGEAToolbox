function [requirerefresh] = callback_SubsampleCells(src, ~, methodoption)

requirerefresh = false;
if nargin < 3, methodoption = []; end

[FigureHandle, sce] = gui.gui_getfigsce(src);

if ~strcmp(gui.myQuestdlg(FigureHandle, 'This function subsamples 50% of cells. Continue?',''), 'Yes'), return; end

if isempty(methodoption)
    answer = gui.myQuestdlg(FigureHandle, 'Select method:', '', ...
        {'Uniform Sampling', ...
        'Geometric Sketching [PMID:31176620]'}, 'Uniform Sampling');
    switch answer
        case 'Uniform Sampling'
            methodoption = 1;
        case 'Geometric Sketching [PMID:31176620]'
            if ~pkg.i_checkpython
                gui.myWarndlg(FigureHandle, 'Python not installed.','');
                return;
            end
            methodoption = 2;
        otherwise
            return;
    end
end

tn = round(sce.NumCells/2);
if methodoption == 1
    rng("shuffle");
    idx = randperm(sce.NumCells);
    ids = idx(1:tn);
elseif methodoption == 2
    if ~gui.gui_showrefinfo('Geometric Sketching [PMID:31176620]', ...
            FigureHandle), return; end
    fw = gui.myWaitbar(FigureHandle);
    try
        Xn = log1p(sc_norm(sce.X))';
        if issparse(Xn)
            Xn = full(Xn);
        end
        [~, Xn] = pca(Xn, 'NumComponents', 300);
        ids = run.py_geosketch(Xn, tn);
    catch ME
        gui.myWaitbar(FigureHandle, fw);
        gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
        return;
    end
    gui.myWaitbar(FigureHandle, fw);
end

if ~isempty(ids)
    sce = sce.selectcells(ids);
    gui.myGuidata(FigureHandle, sce, src);
    requirerefresh = true;
else
    gui.myErrordlg(FigureHandle, 'Runtime error. No action is taken.');
end
end
