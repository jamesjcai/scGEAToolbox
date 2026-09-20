function i_celltypedatatip(h, sce)
%I_CELLTYPEDATATIP Label a cell scatter's data tips by cell type.
%
%   gui.i_celltypedatatip(h, sce)
%
%   Replaces the default X/Y/Z/Color rows of the main scatter's data tip
%   with the one thing a user wants when they click a cell - what type it
%   is. The embedding coordinates carry no meaning on their own, a UMAP x
%   of 24.03 says nothing, and "Color 3" is the group's index rather than
%   its name.
%
%   The colon is part of the row's label because DATATIPTEXTROW does not
%   add one: it prints the label, a space, then the value.
%
%   An SCE with no annotation reads "Undetermined", which is the answer
%   to the question rather than a failure to answer it. The default tip
%   is kept only for a scatter drawn from a subset of the cells, where
%   the labels do not line up with the points and no row can be trusted.
%
%   Callers with a better label for the points - the cluster labeller,
%   the cell state picker - set DataTipTemplate themselves after the
%   redraw and so overwrite this.
%
% See also GUI.I_APPLYDISPLAY, GUI.I_GSCATTER3.

if nargin < 2 || isempty(h) || isempty(sce), return; end
if ~isscalar(h) || ~all(pkg.i_isvalid(h)), return; end
if ~isprop(h, 'DataTipTemplate'), return; end

n = numel(h.XData);
if n == 0, return; end

celltype = i_celltypelabels(sce, n);
if isempty(celltype), return; end

% Cell type names are full of underscores, which TeX would silently turn
% into subscripts, so escape them - but keep TeX. Subtypes are written
% "Fibroblasts_{Inflammatory}" to render as subscripts, and the template
% outlives this call: the cluster labeller and the subtype annotator later
% swap in rows of their own that are TeX-escaped the same way, and an
% Interpreter of 'none' left behind here would show their "_{...}" raw.
celltype = gui.i_escapeunderscore(celltype);
h.DataTipTemplate.DataTipRows = dataTipTextRow('Cell type:', cellstr(celltype));
h.DataTipTemplate.Interpreter = 'tex';
end

function labels = i_celltypelabels(sce, n)
%I_CELLTYPELABELS The cell type of each of N plotted cells, ready to show.
%
%   Empty when the SCE's cells cannot be matched to the N points - a
%   scatter drawn from a subset, say - and nothing can honestly be said
%   about them.
%
%   Numeric input goes through DOUBLE first, because STRING of a sparse
%   array throws and anything summed out of SCE.X is sparse. Text is left
%   to STRING, which takes string arrays, cellstr and categorical as they
%   are - DOUBLE would quietly turn "CD8 T cell" into NaN.
labels = strings(0, 1);

x = sce.c_cell_type_tx;
if isempty(x)
    % Never annotated. Only safe to say so for a plot of the whole SCE.
    if n ~= sce.NumCells, return; end
    labels = repmat("Undetermined", n, 1);
    return;
end
if numel(x) ~= n, return; end

if isnumeric(x) || islogical(x)
    labels = string(full(double(x(:))));
else
    try
        labels = string(x(:));
    catch
        labels = strings(0, 1);
        return;
    end
end
labels = strtrim(labels);
% The constructor's lowercase placeholder is not an answer.
labels(ismissing(labels) | labels == "" | strcmpi(labels, "undetermined")) = "Undetermined";
end
