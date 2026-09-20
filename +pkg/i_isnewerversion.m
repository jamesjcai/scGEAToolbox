function tf = i_isnewerversion(vcandidate, vreference)
%I_ISNEWERVERSION True when vcandidate is a strictly newer version than vreference.
%   TF = I_ISNEWERVERSION(VCANDIDATE, VREFERENCE) compares two dotted
%   version strings such as "26.5.1" or "v26.5.1" component by component.
%   Missing trailing components are treated as zero, so "26.6" is newer
%   than "26.5.9". TF is false when either version cannot be parsed.

tf = false;
a = i_parseversion(vcandidate);
b = i_parseversion(vreference);
if isempty(a) || isempty(b)
    return;
end
n = max(numel(a), numel(b));
a(end+1:n) = 0;
b(end+1:n) = 0;
d = find(a~=b, 1);
tf = ~isempty(d) && a(d)>b(d);
end

function parts = i_parseversion(v)
% Convert a version string such as "26.5.1" into a numeric vector.
parts = [];
if isempty(v)
    return;
end
tok = regexp(char(string(v)), '\d+(\.\d+)*', 'match', 'once');
if isempty(tok)
    return;
end
parts = str2double(strsplit(tok, '.'));
if any(isnan(parts))
    parts = [];
end
end
