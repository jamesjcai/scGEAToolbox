function callback_LASSO(src, ~)

[FigureHandle, sce] = gui.gui_getfigsce(src);

answer = gui.myQuestdlg(FigureHandle, 'Select a dependent variable y. Continue?','');
if ~strcmp(answer,'Yes'), return; end

[thisx, ~] = gui.i_select1state(sce, false, false, false, true, FigureHandle);
if isempty(thisx), return; end
if ~isnumeric(thisx)
    gui.myWarndlg(FigureHandle, 'This function works with continuous variables only.');
    return;
end

[Xt] = gui.i_transformx(sce.X, true, "libsize_log1p", FigureHandle);
if isempty(Xt), return; end

if ~gui.i_resetrngseed(src, [], false), return; end   % cancelled

gui.LassoAnalysisApp(Xt', thisx, sce.g, FigureHandle);

end
