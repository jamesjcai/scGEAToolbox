function [Tm] = i_parsemarkerlist(txt)
%I_PARSEMARKERLIST Parse a marker list text block into a two-column table.
%
%   Tm = pkg.i_parsemarkerlist(txt)
%
%   txt is the text the marker editor returns: one cell type (or subtype) per
%   line, its name and its marker genes separated by a tab,
%
%       T memory cells<TAB>CCR7,SELL,IL7R
%
%   Accepted as well is a line whose first separator is a comma, which is what
%   a pasted CSV row looks like, and a line whose genes are separated by
%   semicolons or spaces rather than commas.
%
%   Tm has one row per non-empty line, Var1 the name and Var2 the marker genes
%   as one upper case comma separated string - the shape both
%   PKG.E_MARKERWEIGHT and PKG.E_DETERMINECELLTYPE expect.
%
%   Lines without a name or without any gene are dropped. A marker list is
%   typed by hand often enough that a stray blank line should not become a
%   nameless cell type that every cluster can be assigned to.
%
%   See also pkg.i_markerweights, pkg.i_readmarkertable,
%   pkg.e_markerweight.

Tm = table(strings(0, 1), strings(0, 1));

if nargin < 1 || isempty(txt), return; end
if iscell(txt) && isscalar(txt) && (iscell(txt{1}) || isstring(txt{1}))
    txt = txt{1};   % what GUI.MYINPUTWIN hands back: {{line; line; ...}}
end

lines = string(txt);
lines = lines(:);
% A single char block with embedded newlines still has to become lines.
lines = splitlines(strjoin(lines, newline));
lines = strtrim(lines);
lines = lines(strlength(lines) > 0);
if isempty(lines), return; end

name = strings(numel(lines), 1);
genes = strings(numel(lines), 1);
for k = 1:numel(lines)
    [name(k), genes(k)] = in_splitline(lines(k));
end

genes = upper(genes);
genes = erase(genes, " ");
genes = replace(genes, [";", "|", sprintf('\t')], ",");
genes = regexprep(genes, ',+', ',');
genes = strip(genes, ',');

keep = strlength(strtrim(name)) > 0 & strlength(genes) > 0;
Tm = table(strtrim(name(keep)), genes(keep));
end

function [name, genes] = in_splitline(line)
% Name first, genes second. The tab wins when there is one, because a name may
% itself contain a comma ("Tumor cells, malignant") while a tab-delimited line
% never means anything else.

name = line;
genes = "";

pos = strfind(line, sprintf('\t'));
if isempty(pos)
    % Space-padded columns are what a tab looks like after a round trip
    % through a text box that expands tabs.
    pos = strfind(line, "  ");
end
if isempty(pos)
    pos = strfind(line, ",");
end
if isempty(pos), return; end

name = strtrim(extractBefore(line, pos(1)));
genes = strtrim(extractAfter(line, pos(1)));
end
