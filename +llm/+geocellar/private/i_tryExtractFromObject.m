function data = i_tryExtractFromObject(text)
data = {};
startPos = strfind(text, "{");
if isempty(startPos); return; end
endPos = i_matchBracket(text, startPos(1), '{', '}');
if endPos < 0; return; end
try
    obj = jsondecode(text(startPos(1):endPos));
catch ME
    fprintf("[hypothesis_agent] jsondecode object failed: %s\n", ME.message);
    return;
end
if ~isstruct(obj); return; end
fields = fieldnames(obj);
for k = 1:numel(fields)
    val = obj.(fields{k});
    if isstruct(val)
        val = num2cell(val);
        if ~isempty(val) && isstruct(val{1}); data = val; return; end
    elseif iscell(val) && ~isempty(val) && isstruct(val{1})
        data = val;
        return;
    end
end
end
