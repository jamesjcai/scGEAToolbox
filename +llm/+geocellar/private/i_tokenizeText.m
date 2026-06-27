function tokens = i_tokenizeText(txt)
txt = lower(char(i_toTextScalar(txt)));
tokens = regexp(txt, "[a-z0-9]+", "match");
if isempty(tokens)
    tokens = {};
end
tokens = unique(string(tokens), "stable");
end
