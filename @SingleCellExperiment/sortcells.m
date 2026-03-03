function obj = sortcells(obj, idx)
obj.X = obj.X(:, idx);
obj.s = obj.s(idx, :);
obj.c = obj.c(idx);
obj = i_applyCellIndex(obj, idx, 'select');
end
