function sc_grnview(A, nodenames, figname, parentfig)
if nargin < 4, parentfig = []; end
if nargin < 3, figname = ''; end
%GRN network viewer
if isa(A, 'digraph') || isa(A, 'graph')
    G = A;
else
    if nargin < 2
        nodenames = string((1:size(A, 1))');
    end
    G = pkg.i_makegraph(A, nodenames);
end
if nargin < 3
    figname = sprintf('nodes (n=%d, red=TF); edges (blue=positive, red=negative)', ...
        G.numnodes);
end
gui.i_singlegraph(G, figname, parentfig);
end
