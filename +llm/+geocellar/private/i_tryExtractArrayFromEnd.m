function data = i_tryExtractArrayFromEnd(text)
data = {};
tok = regexp(text, '\[\s*\{', 'start', 'once');
if isempty(tok); return; end
rbrPositions = fliplr(strfind(text, ']'));
for ii = 1:numel(rbrPositions)
    endCandidate = rbrPositions(ii);
    if endCandidate <= tok; break; end
    data = i_tryDecode(text(tok:endCandidate));
    if ~isempty(data)
        fprintf("[hypothesis_agent] Strategy 4 succeeded (tried %d ] candidates).\n", ii);
        return;
    end
end
end
