function sc_uitabgrpfig_scatter(sce, glist, thisx, xlabelv, parentfig)
%see also: gui.i_violinplot, gui.sc_uitabgrpfig_vioplot

if nargin < 5, parentfig = []; end

[Xt] = gui.i_transformx(sce.X);
if isempty(Xt), return; end
n = length(glist);
y=cell(n,1);
for k=1:n
    y{k} = full(Xt(upper(sce.g) == upper(glist(k)), :));
end
gui.i_scattertabs(y, glist, thisx, xlabelv, parentfig);

