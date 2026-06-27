function txt = i_checkApiError(txt, logFcn, passName)
% Detect API-level error messages returned as response content (not exceptions).
% TAMU/LiteLLM embeds errors in the response body instead of returning HTTP 4xx/5xx.
errPatterns = ["🚫", "litellm.", "Traceback (most recent", "APIConnectionError", ...
               "An unexpected error occurred"];
for k = 1:numel(errPatterns)
    if contains(txt, errPatterns(k))
        preview = extractBefore(txt + " ", min(strlength(txt) + 1, 200));
        logFcn(sprintf("WARNING: %s API error in response: %s", passName, preview));
        txt = "";
        return;
    end
end
end
