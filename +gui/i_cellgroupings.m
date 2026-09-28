function [names, values] = i_cellgroupings(sce, cellIdx)
%I_CELLGROUPINGS Categorical cell annotations that split cells into groups.
%   [NAMES, VALUES] = gui.i_cellgroupings(SCE) returns the groupings with at
%   least two groups: Cell Type, Cluster ID, Batch ID, Cell Cycle Phase and
%   any categorical custom attribute. NAMES is a string column; VALUES{k}
%   is a string column with one label per cell.
%
%   [NAMES, VALUES] = gui.i_cellgroupings(SCE, CELLIDX) counts groups among
%   the cells in logical CELLIDX only, so a grouping that is constant there
%   is left out. VALUES still hold one label per cell of SCE.
%
%   See also gui.i_pickgrncells.

if nargin < 2 || isempty(cellIdx), cellIdx = true(sce.NumCells, 1); end

names = ["Cell Type"; "Cluster ID"; "Batch ID"; "Cell Cycle Phase"];
values = {sce.c_cell_type_tx; sce.c_cluster_id; sce.c_batch_id; ...
    sce.c_cell_cycle_tx};
numBuiltin = numel(names);
attrNames = string(sce.list_cell_attributes(1:2:end));
attrValues = sce.list_cell_attributes(2:2:end);
names = [names; attrNames(:)];
values = [values; attrValues(:)];

% A numeric custom attribute with this many distinct values is a
% measurement, not a grouping. The built-in ones are groupings however
% many groups they have: a clustering can well exceed this.
maxNumericGroups = 50;
keep = false(numel(names), 1);
for k = 1:numel(names)
    v = values{k};
    if numel(v) ~= sce.NumCells, continue; end
    if k > numBuiltin && isnumeric(v) && numel(unique(v)) > maxNumericGroups
        continue;
    end
    v = string(v(:));
    values{k} = v;
    inCells = v(cellIdx);
    keep(k) = numel(unique(inCells)) > 1 && ~all(inCells == "undetermined");
end
names = names(keep);
values = values(keep);
end
