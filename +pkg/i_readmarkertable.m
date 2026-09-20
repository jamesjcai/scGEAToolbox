function [txt, srcfile] = i_readmarkertable(filename)
%I_READMARKERTABLE Read a marker gene file into the editor's text format.
%
%   [txt, srcfile] = pkg.i_readmarkertable(filename)
%
%   Reads a user supplied marker list and returns it as the same text block the
%   marker editor shows - one line per cell type, name and genes separated by a
%   tab - so a file and a paste reach PKG.I_PARSEMARKERLIST by one route.
%
%   Recognized:
%     .xlsx .xls   a sheet whose first column names the type and whose second
%                  column lists its markers. Columns are found by name when
%                  they carry one this toolbox already uses (SubType,
%                  CellType, cellName, Name; PositiveMarkers, Markers,
%                  geneSymbolmore1, Genes, MarkerGenes), otherwise the first
%                  two columns are taken in order.
%     .txt .csv .tsv .text
%                  one line per type: the name, then a tab or the first comma,
%                  then the genes.
%
%   txt is "" when the file holds nothing usable; the caller decides whether
%   that is worth a dialog.
%
%   See also pkg.i_parsemarkerlist, gui.i_getcustomsubtypemarkers.

txt = "";
srcfile = "";
if nargin < 1 || isempty(filename), return; end
filename = char(filename);
if ~exist(filename, 'file')
    error('pkg:i_readmarkertable:FileNotFound', ...
        'Marker file not found: %s', filename);
end
srcfile = string(filename);

[~, ~, ext] = fileparts(filename);
switch lower(ext)
    case {'.xlsx', '.xls', '.xlsm'}
        T = readtable(filename, 'TextType', 'string');
        [name, genes] = in_pickcolumns(T);
    otherwise
        lines = in_readlines(filename);
        Tm = pkg.i_parsemarkerlist(lines);
        if isempty(Tm), return; end
        name = string(Tm.Var1);
        genes = string(Tm.Var2);
end

if isempty(name), return; end

% Several rows for one type - which is how a long marker set is usually
% written down - become one line, in the order they appear in the file.
[uname, ~, back] = unique(name, 'stable');
merged = strings(numel(uname), 1);
for k = 1:numel(uname)
    merged(k) = strjoin(genes(back == k), ",");
end
merged = regexprep(merged, ',+', ',');
merged = strip(merged, ',');

keep = strlength(strtrim(uname)) > 0 & strlength(merged) > 0;
if ~any(keep), return; end

txt = strjoin(uname(keep) + sprintf('\t') + merged(keep), newline);
end

function [name, genes] = in_pickcolumns(T)
% The name column and the marker column, by header when the header is one we
% know and by position when it is not.

name = strings(0, 1);
genes = strings(0, 1);
if isempty(T) || width(T) < 2, return; end

vars = string(T.Properties.VariableNames);
namecol = in_findcol(vars, ["SubType", "CellType", "cellName", "Name", ...
    "celltype", "subtype", "type"]);
genecol = in_findcol(vars, ["PositiveMarkers", "Markers", ...
    "geneSymbolmore1", "Genes", "MarkerGenes", "markers", "genes"]);

if isempty(namecol), namecol = 1; end
if isempty(genecol) || genecol == namecol
    genecol = find((1:width(T)) ~= namecol, 1);
end

name = in_tostring(T{:, namecol});
genes = in_tostring(T{:, genecol});
end

function idx = in_findcol(vars, wanted)
idx = [];
for k = 1:numel(wanted)
    hit = find(strcmpi(vars, wanted(k)), 1);
    if ~isempty(hit), idx = hit; return; end
end
end

function s = in_tostring(col)
s = strtrim(string(col(:)));
s(ismissing(s)) = "";
end

function lines = in_readlines(filename)
% READLINES is R2020b and later; FILEREAD is the fallback that works
% everywhere and costs nothing here, the files being a few kilobytes.

txt = fileread(filename);
lines = splitlines(string(txt));
end
