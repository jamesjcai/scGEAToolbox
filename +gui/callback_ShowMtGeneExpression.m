function callback_ShowMtGeneExpression(src, ~)

FigureHandle = src.Parent.Parent;
sce = guidata(FigureHandle);

idx = startsWith(sce.g, 'mt-', 'IgnoreCase', true);
n = sum(idx);
if n > 0
    [ax, bx] = view(findall(FigureHandle,'type','axes'));
    if n <= 9
        gui.i_markergenespanel(sce.X, sce.g, sce.s, ...
            sce.g(idx), [], 9, ax, bx, 'Mt-genes');
    else
        gui.i_markergenespanel(sce.X, sce.g, sce.s, ...
            sce.g(idx), [], 16, ax, bx, 'Mt-genes');
    end
else
    warndlg('No mt-genes found');
end