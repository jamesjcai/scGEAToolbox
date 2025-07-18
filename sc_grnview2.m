function sc_grnview2(A1, A2, nodenames, figname, parentfig)
    % SC_GRNVIEW2  Display two gene regulatory networks side‑by‑side.
    %
    %   sc_grnview2(A1, A2, nodenames) requires adjacency matrices or graph objects A1, A2 plus node names.
    %   sc_grnview2(..., figname) sets a custom figure title.
    %   sc_grnview2(..., figname, parentfig) embeds viewer in an existing figure.

    % Parse inputs robustly
    p = inputParser;
    addRequired(p, 'A1');
    addRequired(p, 'A2');
    addRequired(p, 'nodenames', @(x) isstring(x) || iscellstr(x));
    addOptional(p, 'figname', '', @ischar);
    addOptional(p, 'parentfig', [], @(x) isempty(x) || isgraphics(x));
    parse(p, A1, A2, nodenames, figname, parentfig);

    % Extract validated inputs
    A1 = p.Results.A1;
    A2 = p.Results.A2;
    nodenames = p.Results.nodenames;
    figname = p.Results.figname;
    parentfig = p.Results.parentfig;

    % Ensure nodenames is string vector
    if iscellstr(nodenames)
        nodenames = string(nodenames);
    end
    nodenames = nodenames(:);  % Ensure column vector

    % Construct graph objects
    if isa(A1, 'digraph') || isa(A1, 'graph')
        G1 = A1;
    else
        G1 = pkg.i_makegraph(A1, nodenames);
    end
    if isa(A2, 'digraph') || isa(A2, 'graph')
        G2 = A2;
    else
        G2 = pkg.i_makegraph(A2, nodenames);
    end

    % Set default figure name if none given
    if isempty(figname)
        figname = sprintf('Network comparison (n=%d)', numel(nodenames));
    end

    % Launch the dual-graph viewer GUI
    gui.i_doublegraphs(G1, G2, figname, parentfig);
end

%{
function sc_grnview2(A1, A2, nodenames, figname, parentfig)
%GRN network viewer (two networks)
if nargin < 5, parentfig = []; end
if nargin < 4, figname = ''; end
if nargin < 3, error('USAGE: sc_grnview2(A1,A2,g)'); end
G1 = pkg.i_makegraph(A1, nodenames);
G2 = pkg.i_makegraph(A2, nodenames);
gui.i_doublegraphs(G1, G2, figname, parentfig);
end
%}