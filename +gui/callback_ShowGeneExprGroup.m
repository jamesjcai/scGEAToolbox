function callback_ShowGeneExprGroup(src, ~)

    FigureHandle = src.Parent.Parent;
    sce = guidata(FigureHandle);
    [axx, bxx] = view(findall(FigureHandle,'type','axes'));
    [glist] = gui.i_selectngenes(sce, [], FigureHandle);
    if isempty(glist), return; end
    % answer = questdlg("Select the type of expression values","",...
    %     "Raw UMI Counts","Library Size-Normalized",)


    allowunique = false;
    [thisc] = gui.i_select1class(sce, allowunique);
    if isempty(thisc), return; end

    gui.i_feaplotarray(sce, glist, thisc, false, FigureHandle);

end