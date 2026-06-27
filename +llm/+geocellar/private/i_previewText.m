function txt = i_previewText(text, n)
txt = char(i_toTextScalar(text));
if isempty(txt)
    txt = "";
    return;
end
txt = txt(1:min(numel(txt), n));
end
