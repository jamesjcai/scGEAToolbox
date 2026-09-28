function A = sc_grn(X, type, varargin)
%SC_GRN Construct single-cell gene regulatory network (scGRN)
%
%   A = sc_grn(X) uses default type 'pcrnet'.
%   A = sc_grn(X, type) builds the network with the method named TYPE.
%   A = sc_grn(X, type, extra...) passes further arguments to methods
%   that take them ('grnformer' needs tf_idx; 'tn' takes a gene list).
%
%   The methods, with a one-line description and the extra arguments each
%   accepts, are the rows of net.grnmethods(); run it to list them:
%       T = net.grnmethods();
%       disp(T(:, ["Key", "Summary", "ExtraArgs"]))
%
%   X is a genes-by-cells matrix, expected to be already normalised and
%   transformed, as it is on every in-tree path: the GUI runs
%   gui.i_transformx first, and +cli/cmd_grn normalises before calling.
%   No branch normalises for you.
%
%   A is genes-by-genes, indexed in the order of the rows of X. A(i, j) is
%   the edge from gene i to gene j (regulators in rows) for every method;
%   symmetric methods give A(i, j) == A(j, i). net.pcrnet itself returns
%   targets in rows, and its registry rows transpose it.
%
%   See also: net.grnmethods, sc_grnview.

arguments
    X {mustBeNumeric}
    type (1,1) string = "pcrnet"
end
arguments (Repeating)
    varargin
end

catalog = net.grnmethods();
row = find(catalog.Key == lower(type));
if isempty(row)
    error("sc_grn:InvalidType", "Type must be one of: %s", ...
        strjoin(catalog.Key, ", "));
end
if ~isempty(varargin) && catalog.ExtraArgs(row) == ""
    error("sc_grn:UnexpectedArguments", ...
        "Method '%s' takes no arguments after the type; got %d.", ...
        catalog.Key(row), numel(varargin));
end

A = catalog.Build{row}(X, varargin);
end
