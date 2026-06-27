function msg = i_buildUserMessage(topic, retrievedText, seed, past)
% Truncate retrieved text to avoid gateway timeouts from very long model
% thinking sessions on TAMU and similar proxies with short SSE timeouts.
maxChars = 8000;
if strlength(retrievedText) > maxChars
    retrievedText = extractBefore(retrievedText, maxChars) + ...
        sprintf("\n... [truncated to %d chars]", maxChars);
end
msg = sprintf("Research topic: %s\n\nRetrieved GEO datasets:\n%s\n", topic, retrievedText);
if ~isempty(past)
    msg = msg + sprintf("\nPreviously analysed (avoid repeating):\n%s\n", past);
end
if ~isempty(seed)
    msg = msg + sprintf("\nSeed hypothesis from user:\n%s\n", seed);
end
msg = msg + newline + ...
    "Return ONLY a valid JSON array of hypothesis objects that follows the required schema exactly.";
end
