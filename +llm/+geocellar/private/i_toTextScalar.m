function txt = i_toTextScalar(val)
if nargin == 0 || isempty(val)
    txt = "";
    return;
end
if isstring(val)
    txt = strjoin(val(:)', newline);
elseif ischar(val)
    txt = string(val);
elseif iscell(val)
    parts = cell(1, 0);
    for k = 1:numel(val)
        part = char(i_toTextScalar(val{k}));
        if ~isempty(strtrim(part))
            parts{end+1} = part; %#ok<AGROW>
        end
    end
    txt = string(strjoin(parts, newline));
elseif isnumeric(val) || islogical(val)
    txt = strjoin(string(val(:))', ", ");
else
    try
        txt = string(val);
    catch
        % fall back to empty if conversion to string is unsupported
        txt = "";
    end
end
txt = string(txt);
if numel(txt) > 1
    txt = strjoin(txt(:)', newline);
end
end
