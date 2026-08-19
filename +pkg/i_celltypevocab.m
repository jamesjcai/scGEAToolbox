function vocab = i_celltypevocab()
% I_CELLTYPEVOCAB - Normalized vocabulary of known cell type names.
%
%   vocab = pkg.i_celltypevocab()
%
% Returns a string array of cell type and cell subtype names collected from the
% bundled PanglaoDB tables, normalized by pkg.i_normalizetypename. The result is
% cached in a persistent variable because the spreadsheets are read from disk.
%
% Used by pkg.i_guesscelltypecol to recognize a metadata column whose values are
% cell type labels, regardless of what the column is named.
%
% see also: pkg.i_guesscelltypecol, pkg.i_normalizetypename

persistent cachedvocab

if ~isempty(cachedvocab)
    vocab = cachedvocab;
    return;
end

pw1 = fileparts(mfilename('fullpath'));
sources = {'celltypes.xlsx', 'CellType'; 'cellsubtypes.xlsx', 'SubType'};
collected = cell(size(sources, 1), 1);

for k = 1:size(sources, 1)
    collected{k} = in_readnames( ...
        in_locate(pw1, sources{k, 1}), sources{k, 2});
end

vocab = unique(pkg.i_normalizetypename(vertcat(collected{:})));
vocab(strlength(vocab) < 2) = [];
cachedvocab = vocab;
end

function pth = in_locate(pw1, filenm)
% Resolve an asset both in the source tree and in a deployed layout.
pth = fullfile(pw1, '..', 'assets', 'PanglaoDB', filenm);
if ~exist(pth, 'file')
    pth = which(filenm);
end
end

function names = in_readnames(pth, colname)
names = strings(0, 1);
if isempty(pth) || ~exist(pth, 'file'), return; end
try
    sheets = sheetnames(pth);
catch
    % Unreadable spreadsheet; the other source may still provide names.
    return;
end
persheet = cell(numel(sheets), 1);
for k = 1:numel(sheets)
    persheet{k} = strings(0, 1);
    try
        t = readtable(pth, 'FileType', 'spreadsheet', 'Sheet', sheets(k));
    catch
        continue; % skip sheets that are not marker tables
    end
    if ismember(colname, t.Properties.VariableNames)
        persheet{k} = string(t.(colname));
    end
end
names = vertcat(persheet{:});
end
