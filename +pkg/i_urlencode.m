function urlOut = i_urlencode(urlIn)
%I_URLENCODE Percent-encode text for a URL query, as URLENCODE did.
%   urlOut = pkg.i_urlencode(urlIn)
%
%   Letters, digits and - _ . * are kept, a space becomes '+', and every
%   other character becomes %XX for each of its UTF-8 bytes. MATLAB's own
%   URLENCODE is unsupported and excluded from the MATLAB Runtime, so a
%   compiled scgeatool cannot call it.
%
%   See also pkg.i_urldecoding.

arguments
    urlIn {mustBeText}
end

urlIn = char(urlIn);
if isempty(urlIn)
    urlOut = '';
    return;
end
bytes = unicode2native(urlIn, 'UTF-8');
safe = (bytes >= 'a' & bytes <= 'z') | (bytes >= 'A' & bytes <= 'Z') | ...
    (bytes >= '0' & bytes <= '9') | ismember(bytes, uint8('-_.*'));
parts = cell(1, numel(bytes));
parts(safe) = cellstr(char(bytes(safe))');
parts(bytes == ' ') = {'+'};
other = ~safe & bytes ~= ' ';
parts(other) = cellstr(strcat('%', dec2hex(bytes(other), 2)))';
urlOut = [parts{:}];
end
