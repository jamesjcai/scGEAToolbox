function data = i_tryExtractArray(text)
data = {};
tok = regexp(text, '\[\s*\{', 'start', 'once');
if isempty(tok)
    % fprintf("[hypothesis_agent] i_tryExtractArray: no [{ pattern found\n");
    return;
end
startPos = tok;
endPos = i_matchBracket(text, startPos, '[', ']');
if endPos < 0
    % fprintf("[hypothesis_agent] i_tryExtractArray: bracket matching failed at pos %d\n", startPos);
    return;
end
try
    parsed = jsondecode(text(startPos:endPos));
    if isstruct(parsed); parsed = num2cell(parsed); end
    if iscell(parsed) && ~isempty(parsed)
        data = parsed;
    end
catch ME
    % fprintf("[hypothesis_agent] jsondecode array failed: %s\n", ME.message);
end
end
