function [bestname, bestscore, ranked] = i_guesscelltypecol(t, ncells)
% I_GUESSCELLTYPECOL - Rank metadata columns by how likely they hold cell types.
%
%   [bestname, bestscore, ranked] = pkg.i_guesscelltypecol(t)
%   [bestname, bestscore, ranked] = pkg.i_guesscelltypecol(t, ncells)
%
% Inputs:
%   t      - per-cell metadata table, one row per cell (e.g. Seurat meta.data)
%   ncells - expected number of cells (default: height(t))
%
% Outputs:
%   bestname  - name of the highest scoring column, "" if none qualifies
%   bestscore - score of that column, 0 to 100; 70 or above is confident
%   ranked    - table of all qualifying columns sorted by Score descending, with
%               variables Name, Score, NumTypes and Examples
%
% The score combines two independent signals. The first is the column name,
% matched against common cell type field names, with QC, batch and sample fields
% blacklisted. The second is the column content: the fraction of distinct labels
% found in the PanglaoDB cell type vocabulary. Because the content term does not
% depend on the name, a column called "foo" holding "B cells" and "Monocytes" is
% still recognized, which a name lookup table can never do.
%
% Example:
%   [~, metadata] = sc_readrdsfile("pbmc.rds");
%   [colname, score] = pkg.i_guesscelltypecol(metadata)
%
% see also: gui.i_pickcelltypecolumn, pkg.i_celltypevocab, sc_readrdsfile

if nargin < 2 || isempty(ncells), ncells = height(t); end

bestname = "";
bestscore = 0;
ranked = table(strings(0, 1), zeros(0, 1), zeros(0, 1), strings(0, 1), ...
    'VariableNames', {'Name', 'Score', 'NumTypes', 'Examples'});

if ~istable(t) || isempty(t) || height(t) ~= ncells, return; end

vocab = pkg.i_celltypevocab();
vocablong = vocab(strlength(vocab) >= 4);

% A cell type column has many fewer distinct values than cells. The cap keeps
% barcode-like and free-text columns out without excluding finely annotated
% atlases, which rarely exceed a couple hundred labels.
maxtypes = max(2, min(200, ceil(ncells/5)));

varnames = string(t.Properties.VariableNames).';
nvar = numel(varnames);
isok = false(nvar, 1);
score = zeros(nvar, 1);
numtypes = zeros(nvar, 1);
examples = strings(nvar, 1);

for k = 1:nvar
    v = t.(varnames(k));
    if ~(isstring(v) || iscellstr(v) || iscategorical(v)), continue; end
    if size(v, 2) ~= 1, continue; end

    u = string(v);
    u = u(~ismissing(u));
    u = strtrim(u);
    u = unique(u(strlength(u) > 0));
    nu = numel(u);
    if nu < 2 || nu > maxtypes, continue; end

    % All-numeric labels are cluster ids stored as text, not cell types.
    if all(~isnan(str2double(u))), continue; end

    isok(k) = true;
    numtypes(k) = nu;
    examples(k) = strjoin(u(1:min(4, nu)), ", ");
    score(k) = in_combinescore(in_namescore(varnames(k)), ...
        in_vocabfrac(u, vocab, vocablong));
end

if ~any(isok), return; end

ranked = table(varnames(isok), score(isok), numtypes(isok), examples(isok), ...
    'VariableNames', {'Name', 'Score', 'NumTypes', 'Examples'});
ranked = sortrows(ranked, 'Score', 'descend');
bestname = ranked.Name(1);
bestscore = ranked.Score(1);
end

function sc = in_combinescore(namescore, vocabfrac)
% Combine the two signals so that either one alone can carry a column, with a
% bonus when they agree. Averaging them instead would sink a decisively named
% column such as predicted.celltype.l2, whose abbreviated labels ("CD14 Mono",
% "pDC") are absent from the vocabulary. A blacklisted name vetoes instead: the
% column can still be offered for manual selection but never taken on its own.

vocabscore = 100*vocabfrac;
if namescore < 0
    sc = max(0, vocabscore + namescore);
    return;
end
sc = min(100, max(namescore, vocabscore) + 0.25*min(namescore, vocabscore));
end

function sc = in_namescore(rawname)
% Score a column name from -50 (known to be something else) to 100.
nm = lower(regexprep(char(rawname), '[^A-Za-z0-9]', ''));

exclude = {'origident', 'barcode', 'cellid', 'cellbarcode', 'cellname', ...
    'sample', 'sampleid', 'samplename', 'patient', 'patientid', 'donor', ...
    'donorid', 'subject', 'subjectid', 'batch', 'batchid', 'library', ...
    'condition', 'treatment', 'group', 'timepoint', 'time', 'stage', ...
    'tissue', 'organ', 'region', 'sex', 'gender', 'age', 'disease', ...
    'status', 'phase', 'sscore', 'g2mscore', 'doublet', 'doubletscore', ...
    'scrubletscore', 'predicteddoublet', 'seuratclusters'};
excludeprefix = {'ncount', 'nfeature', 'ngene', 'nmolecules', 'percent', 'pct', 'rna'};

if any(strcmp(nm, exclude)) || startsWith(nm, excludeprefix)
    sc = -50;
elseif any(strcmp(nm, {'celltype', 'celltypes', 'ctype', 'ctypes'}))
    sc = 100;
elseif contains(nm, 'celltype')
    sc = 90; % predictedcelltypel2, majorcelltype, finalcelltype, ...
elseif any(strcmp(nm, {'annotation', 'annotations', 'annot', 'anno'}))
    sc = 80;
elseif contains(nm, 'annotation') || contains(nm, 'celltypist')
    sc = 75; % seuratannotations, cellannotation, majorannotation, ...
elseif contains(nm, 'subtype') || contains(nm, 'lineage') || ...
        contains(nm, 'celllabel') || contains(nm, 'clusterlabel') || ...
        contains(nm, 'clustername')
    sc = 65;
elseif any(strcmp(nm, {'ident', 'idents', 'identity', 'activeident'}))
    sc = 55;
elseif contains(nm, 'label') || contains(nm, 'class') || contains(nm, 'type')
    sc = 45;
else
    sc = 0;
end
end

function frac = in_vocabfrac(u, vocab, vocablong)
% Fraction of the distinct labels in u that name a known cell type.
frac = 0;
if isempty(vocab), return; end

un = pkg.i_normalizetypename(u);
hit = false(numel(un), 1);
for k = 1:numel(un)
    x = un(k);
    if strlength(x) < 2, continue; end
    if any(vocab == x)
        hit(k) = true;
    elseif strlength(x) >= 4
        % Substring match both ways so that "cd8 t cells" matches "t cells"
        % and "t cell" matches "t cells". Short tokens are excluded from this
        % test because they match almost anything.
        hit(k) = any(contains(vocablong, x)) || any(contains(x, vocablong));
    end
end
frac = mean(hit);
end
