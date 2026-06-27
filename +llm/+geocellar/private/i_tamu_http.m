function [txt, stop] = i_tamu_http(chat, msgs, cfg, logFcn, passName)
% Direct HTTP call to TAMU API, bypassing openAIChat's SSE parser.
apiKey = char(cfg.OpenAIAPIKey);
url    = [regexprep(char(chat.BaseURL), '/$', '') '/chat/completions'];
model  = char(chat.ModelName);

allMsgs = chat.SystemPrompt;
for k = 1:numel(msgs.Messages)
    m = msgs.Messages{k};
    allMsgs{end+1} = struct('role', char(m.role), 'content', char(m.content)); %#ok<AGROW>
end

reqBody = jsonencode(struct( ...
    'model',      model, ...
    'messages',   {allMsgs}, ...
    'stream',     false, ...
    'max_tokens', 8192 ...
));

opts = weboptions( ...
    'HeaderFields', {'Authorization', ['Bearer ' apiKey]; ...
                     'Content-Type',  'application/json'}, ...
    'RequestMethod', 'post', ...
    'Timeout',       300, ...
    'ContentType',   'text' ...
);

responseText = webwrite(url, reqBody, opts);

if startsWith(strtrim(responseText), '<')
    error('llm:geocellar:tamuHtml', 'TAMU returned HTML error page.');
end

try
    data = jsondecode(responseText);
    txt  = string(data.choices(1).message.content);
    stop = string(data.choices(1).finish_reason);
    logFcn(sprintf("%s [direct]: %d chars.", passName, strlength(txt)));
    return;
catch
    % response was streaming/SSE, not a single JSON; fall through to SSE parser below
end

txt  = i_parse_sse(responseText);
stop = "stop";
logFcn(sprintf("%s [direct SSE]: %d chars.", passName, strlength(txt)));
end
