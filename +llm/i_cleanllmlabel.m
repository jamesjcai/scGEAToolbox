function lbl = i_cleanllmlabel(raw, maxlen)
%I_CLEANLLMLABEL Reduce a chat model's reply to a single short label.
%
%   lbl = LLM.I_CLEANLLMLABEL(raw)
%   lbl = LLM.I_CLEANLLMLABEL(raw, maxlen)
%
%   Local models pad short answers with reasoning traces, chat-template
%   markers, markdown and preambles even when told not to. This strips those
%   and keeps the first non-empty line.
%
%   A reply longer than MAXLEN characters (default 60) is returned as
%   "Unknown" rather than truncated. A model that ignored "name only" and
%   wrote a sentence has not given a label, and half a sentence is a worse
%   answer than an honest refusal.
%
%   Returns "" only when the reply had no content at all, which the caller
%   should treat as a failed attempt rather than as a label.
%
%   See also LLM.E_CELLTYPEANNO, LLM.I_VOTELABELS.

arguments
    raw
    maxlen (1,1) double {mustBeInteger, mustBePositive} = 60
end

txt = string(raw);
if ismissing(txt) || strlength(txt) == 0
    lbl = "";
    return
end

% Reasoning traces first.
txt = regexprep(txt, '<think>[\s\S]*?</think>', '');

% Then the chat-template role HEADER, taken together with the role word it
% introduces. Removing only the <|...|> marker would leave "assistant" alone
% on the first line, and "first non-empty line" would then return that as the
% cell type - a plausible-looking label that is purely an artifact.
txt = regexprep(txt, ...
    '^\s*(<\|im_start\|>|<\|start_header_id\|>)\s*\w+(\s*<\|[^|]*\|>)?\s*', '');
txt = regexprep(txt, '<\|[^|]*\|>', '');
txt = erase(txt, ["**", "`", "#"]);

lines = strtrim(splitlines(txt));
lines = lines(strlength(lines) > 0);
if isempty(lines)
    lbl = "";
    return
end

lbl = lines(1);
lbl = regexprep(lbl, '^\s*(answer|cell type|celltype)\s*[:\-]\s*', '', 'ignorecase');
lbl = strtrim(erase(lbl, [char(34), "."]));

if strlength(lbl) > maxlen
    lbl = "Unknown";
end
end
