function callback_CompareGeneNetwork(src, ~)
[FigureHandle, sce] = gui.gui_getfigsce(src);

[i1, i2, cL1, cL2] = gui.i_select2smplgrps(sce, false, FigureHandle);
if isscalar(i1) || isscalar(i2)
    if i1 == 0 || i2 == 0, return; end
end

[glist] = gui.i_selectngenes(sce,[],FigureHandle);
if isempty(glist), return; end

[y, i] = ismember(glist, sce.g);
if ~all(y), error('Selected gene(s) not in the gene list of data.'); end
fprintf("% s\n", glist)

methods = {'PCR (PC Regression)', ...
           'Chatterjee Xi Correlation', ...
           'Pearson Correlation', ...
           'Distance Correlation', ...
           'Mutual Information', ...
           'GENIE3 (Random Forest)'};
methodkeys = {'pcrnet', 'xicor', 'pearson', 'distcorr', 'mi', 'genie3'};

[sel, ok] = gui.myListdlg(FigureHandle, methods, ...
    'Select GRN construction method', 1, false);
if ~ok, return; end
methodkey = methodkeys{sel};

[Xt] = gui.i_transformx(sce.X, true, 5, FigureHandle);
if isempty(Xt), return; end

x1 = Xt(i, i1);
x2 = Xt(i, i2);

fw = gui.myWaitbar(FigureHandle);
try
    A1 = sc_grn(x1, methodkey);
    A2 = sc_grn(x2, methodkey);
catch ME
    gui.myWaitbar(FigureHandle, fw);
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    return;
end
gui.myWaitbar(FigureHandle, fw);

drawnow;
stitle = sprintf('%s vs. %s', cL1{1}, cL2{1});
try
    sc_grnview2(A1, A2, glist, stitle, FigureHandle);
catch ME
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
end
end
