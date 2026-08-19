function callback_ShowGeneExpr(src, ~)

[FigureHandle, sce] = gui.gui_getfigsce(src);

% sce2 = pkg.upgradeSCE(sce)

[axx, bxx] = view(findall(FigureHandle,'type','axes'));
[glist] = gui.i_selectngenes(sce, [], FigureHandle);
if isempty(glist)
    % gui.myHelpdlg(FigureHandle, 'No genes is selected or found.');
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
    gui.sc_uitabgrpfig_expplot(y, glist, sce.s, FigureHandle, [axx, bxx]);
    gui.myWaitbar(FigureHandle, fw);
end
