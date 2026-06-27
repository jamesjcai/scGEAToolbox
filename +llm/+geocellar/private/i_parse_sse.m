function txt = i_parse_sse(sseText)
txt   = "";
lines = strsplit(char(sseText), newline);
for k = 1:numel(lines)
    line = strtrim(lines{k});
    if ~startsWith(line, 'data: '); continue; end
    payload = line(7:end);
    if strcmp(payload, '[DONE]'); continue; end
    try
        ev    = jsondecode(payload);
        delta = ev.choices(1).delta;
        if isfield(delta, 'content') && ~isempty(delta.content)
            txt = txt + string(delta.content);
        end
    catch
        % Skip malformed or HTML-injected events
    end
end
end
