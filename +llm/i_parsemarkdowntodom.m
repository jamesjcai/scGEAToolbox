function domBlocks = i_parsemarkdowntodom(text)
% I_PARSEMARKDOWNTODOM  Convert a markdown string to mlreportgen DOM blocks.
%
% Supported syntax:
%   Headings       : # H1  through  ###### H6
%   Bold+Italic    : ***text***
%   Bold           : **text**
%   Italic         : *text*
%   Ordered list   : 1. item  2. item  ...
%   Unordered list : - item   or  * item   or  + item
%   Paragraph      : any other non-empty line
%
% Returns a cell array of mlreportgen.dom objects ready to append to a doc:
%   domBlocks = llm.i_parsemarkdowntodom(mdText);
%   for k = 1:numel(domBlocks)
%       append(doc, domBlocks{k});
%   end

import mlreportgen.dom.*

domBlocks = {};
if isempty(text), return; end

% Normalise input to char
if isstring(text), text = char(text); end

lines = splitlines(text);
% Ensure cell array (splitlines on char returns cell, but guard anyway)
if ischar(lines), lines = {lines}; end

n = numel(lines);
i = 1;
while i <= n
    line = strtrim(lines{i});
    if isempty(line), i = i + 1; continue; end

    % ---- Heading  (#, ##, …, ######) -----------------------------------
    hm = regexp(line, '^(#{1,6})\s+(.*)$', 'tokens', 'once');
    if ~isempty(hm)
        level = min(numel(hm{1}), 6);
        h = Heading(level);
        i_appendinline(h, hm{2});
        domBlocks{end+1} = h; %#ok<AGROW>
        i = i + 1;
        continue;
    end

    % ---- Ordered list  (1. item) ----------------------------------------
    if ~isempty(regexp(line, '^\d+\.\s+', 'once'))
        ol = OrderedList();
        while i <= n
            cur = strtrim(lines{i});
            om = regexp(cur, '^\d+\.\s+(.*)$', 'tokens', 'once');
            if isempty(om), break; end
            li = ListItem();
            i_appendinline(li, om{1});
            append(ol, li);
            i = i + 1;
        end
        domBlocks{end+1} = ol; %#ok<AGROW>
        continue;
    end

    % ---- Unordered list  (- item  /  * item  /  + item) -----------------
    if ~isempty(regexp(line, '^[-*+]\s+', 'once'))
        ul = UnorderedList();
        while i <= n
            cur = strtrim(lines{i});
            um = regexp(cur, '^[-*+]\s+(.*)$', 'tokens', 'once');
            if isempty(um), break; end
            li = ListItem();
            i_appendinline(li, um{1});
            append(ul, li);
            i = i + 1;
        end
        domBlocks{end+1} = ul; %#ok<AGROW>
        continue;
    end

    % ---- Regular paragraph ----------------------------------------------
    p = Paragraph();
    i_appendinline(p, line);
    domBlocks{end+1} = p; %#ok<AGROW>
    i = i + 1;
end
end

% -------------------------------------------------------------------------
function i_appendinline(container, text)
% Append inline-formatted Text objects (bold / italic / plain) to a
% DOM container (Heading, Paragraph, or ListItem).
import mlreportgen.dom.*
segs = i_parseinline(text);
if ~iscell(segs), segs = {segs}; end
for k = 1:numel(segs)
    append(container, segs{k});
end
end

% -------------------------------------------------------------------------
function content = i_parseinline(text)
% Split a text string into a mix of plain and formatted Text DOM objects.
% Returns a single Text object or a cell array of Text objects.
import mlreportgen.dom.*

if isempty(text)
    content = Text('');
    return;
end

% Match ***bold+italic***, **bold**, *italic*  (greedy-free, no nested *)
pattern = '(\*\*\*[^*]+\*\*\*|\*\*[^*]+\*\*|\*[^*]+\*)';
[s0, e0, mtch] = regexp(text, pattern, 'start', 'end', 'match');

if isempty(mtch)
    content = Text(text);
    return;
end

content = {};
lastEnd = 0;

for k = 1:numel(mtch)
    % plain text before this match
    if s0(k) > lastEnd + 1
        content{end+1} = Text(text(lastEnd+1 : s0(k)-1)); %#ok<AGROW>
    end
    m = mtch{k};
    if startsWith(m, '***') && endsWith(m, '***')
        t = Text(m(4:end-3)); t.Bold = true; t.Italic = true;
    elseif startsWith(m, '**') && endsWith(m, '**')
        t = Text(m(3:end-2)); t.Bold = true;
    else  % *italic*
        t = Text(m(2:end-1)); t.Italic = true;
    end
    content{end+1} = t; %#ok<AGROW>
    lastEnd = e0(k);
end

% trailing plain text
if lastEnd < numel(text)
    content{end+1} = Text(text(lastEnd+1:end)); %#ok<AGROW>
end

if isscalar(content)
    content = content{1};
end
end
