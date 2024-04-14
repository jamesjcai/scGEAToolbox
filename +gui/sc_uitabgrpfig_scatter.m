function sc_uitabgrpfig_scatter(sce, glist, thisx, xlabelv, parentfig)
%see also: gui.i_violinplot, gui.i_violintabs

if nargin<4, parentfig = []; end

[Xt] = gui.i_transformx(sce.X);
n = length(glist);
y=cell(n,1);
for k=1:n
    y{k} = full(Xt(upper(sce.g) == upper(glist(k)), :));
end
gui.i_scattertabs(y, glist, thisx, xlabelv, parentfig);

