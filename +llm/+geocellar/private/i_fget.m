function val = i_fget(s, fname, default)
if isfield(s, fname) && ~isempty(s.(fname))
    val = s.(fname);
else
    val = default;
end
end
