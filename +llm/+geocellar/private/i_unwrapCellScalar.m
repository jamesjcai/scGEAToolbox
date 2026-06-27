function raw = i_unwrapCellScalar(raw, default)
while iscell(raw)
    if isempty(raw)
        raw = default;
        return;
    end
    raw = raw{1};
end
end
