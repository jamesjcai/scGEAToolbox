function callback_CompareGeneNetwork(src, ~)
FigureHandle = src.Parent.Parent;
sce = guidata(FigureHandle);

[i1, i2] = gui.i_select2grps(sce);
if isscalar(i1) || isscalar(i2)
    if i1 == 0 || i2 == 0, return; end
end

[glist] = gui.i_selectngenes(sce);
if isempty(glist), return; end

[y, i] = ismember(glist, sce.g);
if ~all(y), error('xxx'); end
fprintf("%s\n", glist)
%fw=gui.gui_waitbar;
x1 = sce.X(i, i1);
x2 = sce.X(i, i2);
A1 = sc_pcnet(x1, 3, false, true, true);
A2 = sc_pcnet(x2, 3, false, true, true);
%gui.gui_waitbar(fw);
pause(1)
sc_grnview2(A1, A2, glist);
end
