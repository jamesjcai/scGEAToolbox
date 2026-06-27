function [T] = i_makeattributestable(sce)

n = size(sce.X, 2);   % reliable cell count regardless of c_cell_id storage type

% Normalise c_cell_id → n×1 string array (handles char matrices from old .mat files)
cid = sce.c_cell_id;
if ischar(cid)
    cid = string(cellstr(cid));   % char matrix (n×maxlen) → n strings
else
    cid = string(cid);
end
cid = cid(:);
if numel(cid) ~= n
    cid = string(transpose(1:n));  % fallback: integer labels
end

T = table( ...
    cid, ...
    i_pad(sce.c_batch_id, n), ...
    i_pad(sce.c_cluster_id, n), ...
    i_pad(sce.c_cell_type_tx, n), ...
    i_pad(sce.c_cell_cycle_tx, n), ...
    'VariableNames', {'CellID', 'BatchID', 'ClusterID', 'CellType', 'CellCycle'});
names = sce.list_cell_attributes(1:2:end);
for k = 1:numel(names)
    val = string(sce.list_cell_attributes{2*k});
    if numel(val) == n
        T = addvars(T, val(:), 'NewVariableNames', names{k});
    end
end
end

function col = i_pad(col, n)
if numel(col) ~= n
    col = repmat("", n, 1);
else
    col = col(:);
end
end
