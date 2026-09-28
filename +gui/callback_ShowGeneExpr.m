function callback_ShowGeneExpr(src, ~)

[FigureHandle, sce] = gui.gui_getfigsce(src);

[axx, bxx] = view(findall(FigureHandle,'type','axes'));
[glist] = gui.i_selectngenes(sce, [], FigureHandle);
if isempty(glist)
    return;
end

% answer = gui.myQuestdlg(FigureHandle, "Select the type of expression values","",...
%     "Raw UMI Counts","Library Size-Normalized",)

[Xt] = gui.i_transformx(sce.X,[],[],FigureHandle);
if isempty(Xt), return; end
n = length(glist);
if ~ispref('scgeatoolbox', 'prefcolormapname')
    setpref('scgeatoolbox', 'prefcolormapname', 'autumn');
end

    fw = gui.myWaitbar(FigureHandle);
    y = cell(n, 1);
    for k = 1:n
        y{k} = Xt(sce.g == glist(k), :);
    end    
    try
        gui.sc_uitabgrpfig_expplot(y, glist, sce.s, FigureHandle, [axx, bxx]);
    catch ME
        % The bar is modal on the app; a failed plot used to leave it up.
        gui.myWaitbar(FigureHandle, fw, true);
        gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
        return;
    end
    gui.myWaitbar(FigureHandle, fw);
end
