function txt = i_stripThink(txt)
txt  = char(txt);
inner = regexp(txt, "<think>(.*?)</think>", "tokens", "dotall");
out   = strtrim(regexprep(txt, "<think>.*?</think>", "", "dotall"));
if isempty(out) && ~isempty(inner)
    out = strtrim(char(inner{end}{1}));
end
txt = i_stripChatArtifacts(out);
end
