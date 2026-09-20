function col = groupcol(G, cellType, condition)
%GROUPCOL  Column index of a cell type (and condition) in a glyco-norm struct.
%   COL = GLY.GROUPCOL(G, CELLTYPE) returns the column of G.scores
%   holding the normalized glyco-module scores for CELLTYPE, or 0 when the
%   cell type is not present.
%
%   COL = GLY.GROUPCOL(G, CELLTYPE, CONDITION) selects the column for a
%   specific condition. When G was built without conditions the CONDITION
%   argument is ignored, so callers can pass it unconditionally.
%
%   INPUTS:
%     G         - struct from GLY.NORM
%     cellType  - scalar cell-type label
%     condition - scalar condition label (optional; "" = ignore)
%
%   OUTPUT:
%     col - scalar column index into G.scores, or 0 when there is no match
%
% see also: GLY.NORM, GLY.LRWEIGHT

arguments
    G (1, 1) struct
    cellType (1, 1) string
    condition (1, 1) string = ""
end

if G.hascond && condition ~= ""
    hit = find(G.celltypes == cellType & G.conditions == condition, 1);
else
    hit = find(G.celltypes == cellType, 1);
end

if isempty(hit)
    col = 0;
else
    col = hit;
end

end
