function T = e_cellconditiontable(sce)
%E_CELLCONDITIONTABLE  The per-cell variables cells can be filtered on, in their own types.
%
%   T = pkg.e_cellconditiontable(sce) returns a table with one row per cell:
%   CellID, CellType, Cluster, Batch and CellCycle where SCE has them,
%   the QC metrics LibrarySize, GenesDetected and MitoPercent (the share of
%   counts on genes named "mt-..." in either case, left out when there are
%   none), and then every entry of SCE.LIST_CELL_ATTRIBUTES.
%
%   Unlike PKG.I_MAKEATTRIBUTESTABLE, which turns every column into text
%   for display, numbers stay numbers here, so a condition such as
%   "MitoPercent greater than 10" compares values rather than strings.
%   Categorical and text columns become strings. Columns that are empty, or
%   not one value per cell, are left out.
%
%   See also PKG.E_MATCHCONDITIONS, GUI.CALLBACK_BRUSHCELLSBYCONDITIONS.

arguments
    sce SingleCellExperiment
end

n = sce.NumCells;
names = strings(0, 1);
values = cell(0, 1);

in_add("CellID", sce.c_cell_id);
in_add("CellType", sce.c_cell_type_tx);
in_add("Cluster", sce.c_cluster_id);
in_add("Batch", sce.c_batch_id);
in_add("CellCycle", sce.c_cell_cycle_tx);

% SCE.X is often single; accumulate in double so the percentages are exact.
libSize = full(sum(sce.X, 1, "double")).';
in_add("LibrarySize", libSize);
in_add("GenesDetected", full(sum(sce.X > 0, 1)).');
isMito = startsWith(sce.g, "mt-", IgnoreCase=true);
if any(isMito)
    mitoCounts = full(sum(sce.X(isMito, :), 1, "double")).';
    mitoPercent = 100*mitoCounts./libSize;
    mitoPercent(libSize == 0) = NaN;
    in_add("MitoPercent", mitoPercent);
end

attributes = sce.list_cell_attributes;
for k = 1:2:numel(attributes) - 1
    in_add(string(attributes{k}), attributes{k + 1});
end

names = matlab.lang.makeUniqueStrings(names);
T = table(values{:}, VariableNames=names);

    function in_add(name, value)
        if isempty(value) || numel(value) ~= n
            return;
        end
        value = value(:);
        if isnumeric(value) || islogical(value)
            value = double(full(value));
        else
            value = string(value);
            if all(ismissing(value) | value == "")
                return;
            end
        end
        names(end + 1, 1) = name;
        values{end + 1, 1} = value;
    end
end
