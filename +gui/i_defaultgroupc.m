function [c, cL] = i_defaultgroupc(sce)
%I_DEFAULTGROUPC The grouping the main scatter should be coloured by.
%
%   [c, cL] = gui.i_defaultgroupc(sce)
%
%   SCE.C is what the main scatter colours its points by. A freshly read
%   SCE has it all ones - the constructor's placeholder - so every cell is
%   painted the same colour even when the cells carry a cell type. When
%   SCE.C says nothing, the cell type takes over as the grouping and
%   SCE.C is set to the group indices, so the colorbar, the cluster
%   labels and the group menus all see the grouping the plot is drawn
%   with rather than having to rediscover it.
%
%   An SCE that is already grouped is left alone, with one exception: when
%   SCE.C splits the cells exactly the way the cell type does, the names are
%   used as the labels. SCE.C is a bare index vector and carries no names, so
%   a file saved with SCE.C = findgroups(SCE.C_CELL_TYPE_TX) - which is what
%   every annotation handler in the app writes - used to reopen labelled "1"
%   to "15" instead of by cell type, and on the colorbar in string order at
%   that ("1", "10", "11", ..., "2"). Nothing is lost by preferring the names:
%   the grouping is the same one either way.
%
%   A grouping that is NOT the cell type's - finer clusters, batches, anything
%   else the user put in SCE.C - is still left exactly as it is.
%
%   So is an SCE with no cell type at all, or one whose labels do not line up
%   with its cells. An SCE annotated "undetermined" throughout goes through the
%   same path, which changes no colour - one group either way - but leaves
%   SCE.C and the returned CL consistent.
%
%   SCE is a handle, so SCE.C is updated in place.
%
% See also GUI.I_GSCATTER3, GUI.I_CELLTYPEDATATIP.

c = [];
cL = [];
if nargin < 1 || isempty(sce) || ~isa(sce, 'SingleCellExperiment')
    return;
end
if sce.NumCells == 0
    return;
end

t = sce.c_cell_type_tx;
if ~isempty(t) && numel(t) == sce.NumCells
    % Numeric labels go through DOUBLE first: STRING throws on a
    % sparse array, and anything derived from SCE.X is sparse.
    if isnumeric(t) || islogical(t)
        t = full(double(t(:)));
    end
    t = string(t(:));

    % Either SCE.C says nothing, or it says the same thing the cell type
    % does. In both cases the cell type is the better label source, being
    % the only one of the two that has names.
    if ~i_isgrouping(sce.c, sce.NumCells) || i_samepartition(sce.c, t)
        [c, cL] = findgroups(t);
        sce.c = c;
        return;
    end
end

g = sce.c;
if isnumeric(g) || islogical(g)
    g = full(double(g(:)));
end
[c, cL] = findgroups(string(g));
end

function tf = i_isgrouping(c, n)
%I_ISGROUPING Whether C already splits the N cells into more than one group.
%
%   A placeholder SCE.C - all ones - is not a grouping. Neither is one
%   left empty or out of step with the cells, which cannot colour them.
tf = false;
if isempty(c) || numel(c) ~= n
    return;
end
if isnumeric(c) || islogical(c)
    c = full(double(c(:)));
end
tf = ~isscalar(unique(c));
end

function tf = i_samepartition(c, t)
%I_SAMEPARTITION Whether C and T cut the cells into the same groups.
%
%   Not whether they carry the same values - C is an index vector and T is
%   text - but whether every cell that shares a C value shares a T value and
%   the other way round. That is the case exactly when the number of distinct
%   (C, T) pairs equals the number of distinct C values and the number of
%   distinct T values.
%
%   It has to hold in BOTH directions. A C of 30 subclusters inside 15 cell
%   types gives one cell type per C value, but relabelling by cell type would
%   merge 30 groups into 15 and lose the grouping the file was saved with.

tf = false;
if isnumeric(c) || islogical(c)
    c = full(double(c(:)));
end
gc = findgroups(string(c(:)));
gt = findgroups(t(:));
if isempty(gc) || isempty(gt), return; end

npair = size(unique([gc(:), gt(:)], 'rows'), 1);
tf = npair == max(gc) && npair == max(gt);
end
