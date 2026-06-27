function val = i_str(s, fname, default)
if nargin < 3; default = ""; end
raw = i_fget(s, fname, default);
if iscell(raw)
    raw = i_unwrapCellScalar(raw, default);
end
if isempty(raw)
    raw = default;
end
if isstring(raw)
    if isempty(raw)
        raw = string(default);
    else
        raw = raw(1);
    end
end
if isnumeric(raw) || islogical(raw)
    if isempty(raw)
        raw = string(default);
    else
        raw = string(raw(1));
    end
end
if isstruct(raw)
    raw = string(default);
end
val = char(raw);
end
