function sce = i_addh5adobsattribs(sce, filenm)
%I_ADDH5ADOBSATTRIBS  Keep an .h5ad file's remaining obs columns on the SCE.
%
%   sce = pkg.i_addh5adobsattribs(sce, filenm)
%
% SC_READH5ADFILE takes three things out of the obs table -- the barcodes,
% a batch column and a cell type column -- and every other column was
% dropped on import. Those are the per-cell annotations the file came with:
% ontology terms, sample accessions, QC counts, whatever the depositor
% recorded. They are kept here as named cell attributes.
%
% A column is skipped when it does not have one value per cell, when the
% SCE already carries the same values as its cell IDs, batch or cell type,
% or when an attribute of that name is already set. The check against what
% is already on the object is by value rather than by name, so it holds
% whichever of its aliases SC_READH5ADFILE matched.
%
% Categorical columns then show up as grouping variables throughout the
% GUI, and numeric ones as per-cell measurements; PKG.I_ISGROUPINGVAR is
% what tells the two apart.
%
% See also SC_READH5ADFILE, PKG.I_READH5ADOBS, PKG.I_ISGROUPINGVAR

arguments
    sce (1,1) SingleCellExperiment
    filenm (1,1) string
end

[names, values] = pkg.i_readh5adobs(filenm);
if isempty(names), return; end

% The batch and cell type columns are already on the object, under their
% own fields. Which ones those are comes from PKG.I_H5ADOBSROLES, the same
% answer SC_READH5ADFILE used to read them, rather than from comparing
% values: a caller that replaces unlabeled cells with "undetermined" has
% changed them, so the cell type column would no longer look equal to
% c_cell_type_tx and would be kept a second time under its own name.
roles = pkg.i_h5adobsroles(names);
consumed = [roles.batch, roles.celltype];
consumed = consumed(consumed > 0);

for k = 1:numel(names)
    if ismember(k, consumed), continue; end
    v = values{k};
    if numel(v) ~= sce.NumCells, continue; end
    if sce.hasCellAttribute(char(names(k))), continue; end
    % The barcodes are the DataFrame index and already the cell IDs, but a
    % file can also repeat them as an ordinary column under another name.
    if ~isempty(sce.c_cell_id) && numel(sce.c_cell_id) == numel(v) && ...
            isequal(string(v(:)), string(sce.c_cell_id(:)))
        continue;
    end
    sce.setCellAttribute(char(names(k)), v);
end
end
