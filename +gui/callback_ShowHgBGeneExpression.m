function callback_ShowHgBGeneExpression(src,~)

FigureHandle=src.Parent.Parent;
sce=guidata(FigureHandle);

idx1 = startsWith(sce.g, 'Hba-', 'IgnoreCase', true);
idx2 = startsWith(sce.g, 'Hbb-', 'IgnoreCase', true);
idx3= strcmpi(sce.g,"Alas2");
idx=idx1|idx2|idx3;

if any(idx)
    ttxt = sprintf("%s+", sce.g(idx));
    ci = sum(sce.X(idx, :), 1);
    figure;
    gui.i_stemscatter(sce.s,ci);
    title(ttxt);
else
    warndlg('No HgB-genes found');
end
