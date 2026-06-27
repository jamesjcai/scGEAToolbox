function [txt, stop] = i_generate(chat, msgs, logFcn, passName, cfg)
if nargin < 5; cfg = []; end
txt  = "";
stop = "";

isTAMU = ~isempty(cfg) && isstruct(cfg) && ...
    contains(lower(string(cfg.OpenAIBaseURL)), "tamu");

% For TAMU, use direct non-streaming HTTP first. The SSE/StreamFun path
% has a ~30 s gateway timeout that fires mid-think; stream=false requests
% are given the full weboptions Timeout (300 s) instead.
if isTAMU
    try
        [txt, stop] = i_tamu_http(chat, msgs, cfg, logFcn, passName);
        txt = string(i_stripThinkInline(char(txt)));
        txt = i_checkApiError(txt, logFcn, passName);
        logFcn(sprintf("%s: received %d chars.", passName, strlength(txt)));
        return;
    catch ME
        logFcn(sprintf("%s direct HTTP failed (%s) — falling back to openAIChat.", ...
            passName, ME.message));
    end
end

try
    [txt, stop] = generate(chat, msgs);
    txt = string(i_stripThinkInline(char(txt)));
    txt = i_checkApiError(txt, logFcn, passName);
    logFcn(sprintf("%s: received %d chars.", passName, strlength(txt)));
    if strlength(txt) == 0
        logFcn(sprintf("WARNING: %s returned empty response. Check API key and base URL.", passName));
    end
catch ME
    logFcn(sprintf("ERROR in %s: %s", passName, ME.message));
    rethrow(ME);
end
end
