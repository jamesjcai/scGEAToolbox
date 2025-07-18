function sc_grnview(A, nodenames, figname, parentfig)
    % SC_GRNVIEW  Display a gene regulatory network as a graph GUI.
    %
    %   sc_grnview(A) displays graph from adjacency or graph object A.
    %   sc_grnview(A, nodenames) attaches names to numeric nodes.
    %   sc_grnview(A, nodenames, figname) uses custom figure title.
    %   sc_grnview(A, nodenames, figname, parentfig) embeds in parent GUI.

    % Parse inputs robustly
    p = inputParser;
    addRequired(p, 'A');
    addOptional(p, 'nodenames', "", @(x) isstring(x) || iscellstr(x) || isempty(x));
    addOptional(p, 'figname', '', @ischar);
    addOptional(p, 'parentfig', [], @(x) isempty(x) || isgraphics(x));
    parse(p, A, nodenames, figname, parentfig);

    A = p.Results.A;
    nodenames = p.Results.nodenames;
    figname = p.Results.figname;
    parentfig = p.Results.parentfig;

    % Build or validate graph
    if isa(A, 'digraph') || isa(A, 'graph')
        G = A;
    else
        if isempty(nodenames)
            nodenames = string(1:size(A, 1))';
        elseif iscellstr(nodenames)
            nodenames = string(nodenames);
        end
        G = pkg.i_makegraph(A, nodenames);
    end

    % Default figure title
    if isempty(figname)
        figname = sprintf('Nodes (n=%d, red=TF); edges (blue=positive, red=negative)', ...
                          numnodes(G));
    end

    % Launch GUI viewer
    gui.i_singlegraph(G, figname, parentfig);
end

%{
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
%}
