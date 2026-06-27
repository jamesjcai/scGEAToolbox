function txt = i_stripChatArtifacts(txt)
txt = char(txt);

% Strip leading chat-template role headers that some APIs/models include in
% their response body (e.g. kimi-k2.5 on NVIDIA: "<|im_start|>assistant\n...").
% Without this, the stop-marker logic below would cut at position 1, yielding "".
[s, e] = regexp(txt, '^(<\|im_start\|>|<\|start_header_id\|>)\s*\w+(\s*<\|[^|]*\|>)?\s*', ...
    'start', 'end', 'once');
if ~isempty(s)
    txt = strtrim(txt(e+1:end));
end

% Drop everything after common chat-template stop markers emitted by some
% local models (notably Qwen-family Ollama models).
markers = ["<|endoftext|>", "<|im_start|>", "<|start_header_id|>", "</s>", "<|eot_id|>"];
cutPos = inf;
for j = 1:numel(markers)
    idx = strfind(txt, char(markers(j)));
    if ~isempty(idx)
        cutPos = min(cutPos, idx(1));
    end
end
if isfinite(cutPos)
    txt = txt(1:cutPos-1);
end

% Remove any leftover inline special tokens.
txt = regexprep(txt, "<\|[^>]*\|>", "");

% Some local models begin generating a new chat turn inline.
chatTurnIdx = regexp(txt, '(\r?\n|\s)(user|assistant|system)\s*[:>]', 'once');
if ~isempty(chatTurnIdx)
    txt = txt(1:chatTurnIdx-1);
end

txt = strtrim(txt);
end
