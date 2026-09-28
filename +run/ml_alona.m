function [T] = ml_alona(X, genelist, clusterid, varargin)
%ML_ALONA Score PanglaoDB cell types against each cluster.
%
%   T = ML_ALONA(X, GENELIST, CLUSTERID) returns, per cluster, the ranked
%   cell types and their CTA scores. Markers come from the HUMAN or MOUSE
%   sheet of assets/PanglaoDB/celltypes.xlsx.
%
%   T = ML_ALONA(..., 'species', "mouse") picks the sheet.
%   T = ML_ALONA(..., 'bestonly', false) keeps up to ten rows per cluster
%   instead of the top one.
%
%   This used to be two near-identical functions. ML_ALONA read the text
%   files under external/fun_alona_panglaodb and ML_ALONA_NEW the workbook;
%   the scoring loop was the same in both, and on the same input they
%   ranked the types identically, differing only in the fourth significant
%   figure of the score because the two marker tables are not quite the
%   same vintage. The workbook version survived. Two options went with the
%   text-file version, both already unreachable: 'subtype', whose work
%   SC_CSUBTYPEANNO now does against cellsubtypes.xlsx, and zebrafish,
%   which no caller could request -- CLI.CMD_CELLTYPES rejects any species
%   but human and mouse, and the cell-type callbacks all offer those two.
%   The text files stay where they are; PKG.I_GET_PANGLAODBMARKERS,
%   PKG.E_PRIMARYMARKERS, RUN.R_CLUSTERMOLE and others still read them.
%
% https://alona.panglaodb.se/
% https://academic.oup.com/database/article/doi/10.1093/database/baz046/5427041
% REF: PanglaoDB: a web server for exploration of mouse and human single-cell RNA sequencing data

if isempty(X) || isempty(genelist)
    k = 1;
    T = table("Unknown", 0, 'VariableNames', ...
        {sprintf('C%d_Cell_Type', k), sprintf('C%d_CTA_Score', k)});
    return;
end
if nargin < 3 || isempty(clusterid)
    clusterid = ones(1, size(X, 2));
end
if min(size(clusterid)) ~= 1 || ~isnumeric(clusterid)
    error('CLUSTERID={vector|[]}');
end

p = inputParser;
addRequired(p, 'X', @isnumeric);
addRequired(p, 'genelist', @isstring);
addRequired(p, 'clusterid', @isnumeric);
addOptional(p, 'species', "human", @(x) (isstring(x) | ischar(x)) & ismember(lower(string(x)), ["human", "mouse"]));
addOptional(p, 'bestonly', true, @islogical);
parse(p, X, genelist, clusterid, varargin{:});
species = p.Results.species;
bestonly = p.Results.bestonly;

pth = fullfile(cdgea, 'assets', 'PanglaoDB', 'celltypes.xlsx');

if issparse(X)
    try
        X = full(X);
    catch
        disp('Using sparse input--longer running time is expected.');
    end
end
% warning off
% X=sc_norm(X,"type","deseq");
% warning on
genelist = upper(genelist);

Tm = readtable(pth,'FileType','spreadsheet','Sheet',species);
Tw = pkg.e_markerweight(Tm);

wvalu = Tw.Var2;
wgene = string(upper(Tw.Var1));

[validG, idx1, idx2] = intersect(genelist, wgene);

% Normalize on the full matrix, then subset. Library size has to mean the
% cell's total counts; summing only the matched marker genes makes the
% scale factor depend on which markers happened to intersect genelist.
X = sc_norm(X);
genelist = genelist(idx1);
X = log1p(X(idx1, :));

wvalu = wvalu(idx2);
wgene = wgene(idx2);

celltypev = string(Tm.CellType);
markergenev = string(Tm.PositiveMarkers);
NC = max(clusterid);

S = zeros(length(celltypev), NC);

for j = 1:length(celltypev)
    g = strsplit(markergenev(j), ',');
    g = strtrim(g);
    if strlength(g(end)) == 0
        g = g(1:end-1);
    end
    g = upper(unique(g));
    y = matches(g, genelist);
    if ~any(y), continue; end
    g = g(y);
    Z = zeros(NC, 1);
    ng = zeros(NC, 1);
    for i = 1:length(g)
        gidx = g(i) == genelist;
        wi = wvalu(g(i) == wgene);
        for k = 1:NC
            z = X(gidx, clusterid == k);
            z = mean(z(:));
            Z(k) = Z(k) + z * wi;
            ng(k) = ng(k) + 1;
        end
    end
    for k = 1:NC
        if ng(k) > 0
            S(j, k) = Z(k) ./ nthroot(ng(k), 3);
        else
            S(j, k) = 0;
        end
    end
end
T = table();
for k = 1:NC
    [c, idx] = sort(S(:, k), 'descend');
    if ~isempty(validG)
        T = [T, table(celltypev(idx), c, 'VariableNames', ...
            {sprintf('C%d_Cell_Type', k), sprintf('C%d_CTA_Score', k)})];
    else
        T = [T, table("Unknown", 0, 'VariableNames', ...
            {sprintf('C%d_Cell_Type', k), sprintf('C%d_CTA_Score', k)})];
    end
end
if size(T, 1) > 10
    T = T(1:10, :);
end
if bestonly && size(T, 1) > 1
    T = T(1, :);
end
end
