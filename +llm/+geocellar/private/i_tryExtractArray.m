function data = i_tryExtractArray(text)
data = {};
tok = regexp(text, '\[\s*\{', 'start', 'once');
if isempty(tok)
    return;
end
startPos = tok;
endPos = i_matchBracket(text, startPos, '[', ']');
if endPos < 0
    return;
end
try
    parsed = jsondecode(text(startPos:endPos));
    if isstruct(parsed); parsed = num2cell(parsed); end
    if iscell(parsed) && ~isempty(parsed)
        data = parsed;
    end
catch ME
end
end
