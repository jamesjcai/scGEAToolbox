function callback_TCellExhaustionScores(src, ~)

% https://doi.org/10.1016/j.immuni.2018.04.026
% https://doi.org/10.1038/s41590-022-01262-7

FigureHandle = src.Parent.Parent;
sce = guidata(FigureHandle);
disp('xxx');
%fw=gui.gui_waitbar;
%try
cs = pkg.e_cellscores(sce.X, sce.g, "T Cell Exhaustion");
%catch ME
%    gui.gui_waitbar(fw,true);
%    errordlg(ME.message);
%    return;
%end
%gui.gui_waitbar(fw);
if ~isempty(cs)
    figure;
    gui.i_stemscatter(sce.s, cs);
    zlabel('Score Value')
    title('T Cell Exhaustion Score')
    pause(2)
    if ~istable(cs), cs = table(cs, 'VariableNames', ...
            {'TCellExhaustionScore'}); end
    gui.i_exporttable(cs, false, 'TCellExhaustionScores');
end
end