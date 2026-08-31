function tf = i_hascelltypelabels(sce)
%I_HASCELLTYPELABELS True if SCE has cell-type labels worth preserving.
% The constructor fills c_cell_type_tx with "undetermined", which is not a
% label worth warning about or stashing before an overwrite.

v = string(sce.c_cell_type_tx);
tf = ~isempty(v) && ~all(ismissing(v) | v == "" | v == "undetermined");
end
