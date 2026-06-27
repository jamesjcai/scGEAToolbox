function txt = i_stripThinkInline(txt)
% Strip <think>...</think> blocks (extended thinking) from raw model output.
% Handles partial blocks (timeout mid-think) by falling back to inner text.
txt   = char(txt);
inner = regexp(txt, "<think>(.*?)</think>", "tokens", "dotall");
out   = strtrim(regexprep(txt, "<think>.*?</think>", "", "dotall"));
% Partial <think> block (timeout before </think>): strip it entirely.
out   = strtrim(regexprep(out, "<think>.*$", "", "dotall"));
if isempty(out) && ~isempty(inner)
    out = strtrim(char(inner{end}{1}));
end
txt = out;
end
