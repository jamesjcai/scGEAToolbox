function A = i_causalcccassemble(sce, senders, receivers, o)
%PKG.I_CAUSALCCCASSEMBLE  Shared sender/receiver split, L-R pair and gene
%selection for the causalCCC family of functions.
%
%   A = PKG.I_CAUSALCCCASSEMBLE(sce, senders, receivers, o)
%
%   o is a scalar struct with fields GroupBy, LigandReceptor, SenderGenes,
%   ReceiverGenes, GeneSelection, MaxGenes, MaxPairs, MinFrac, Metadata,
%   MIBins, MIPermutations - the same options PKG.SC_SCE2CAUSALCCC exposes.
%
%   This performs every step both PKG.SC_SCE2CAUSALCCC (which writes the
%   mosaic/state-order/links files) and SC_CAUSALCCCNET (which builds a
%   local approximate network instead) need, so the two agree on which
%   cells, ligand-receptor pairs and intracellular genes are in play. It has
%   no file or figure side effects.
%
%   Not in a "private" folder because both a top-level function
%   (SC_CAUSALCCCNET) and another +pkg function (PKG.SC_SCE2CAUSALCCC) call
%   this - a private folder is visible only to functions in its own
%   immediate parent folder, so it could serve one or the other but not
%   both. SC_CAUSALCCCNET reaches it via the qualified name PKG.I_CAUSALCCCASSEMBLE.
%
%   Returned struct A has fields:
%     isSnd, isRcv   - logical masks into sce's native cell order
%     ord            - [find(isSnd); find(isRcv)], senders-first row order
%     g, gU          - sce.g as string, and its upper-case version
%     Xs, Xr         - full double genes-by-cells matrices for each side
%     pairs          - selected [ligand, receptor] pairs, upper-case
%     ligands, receptors - unique ligands/receptors among pairs, 'stable'
%     sndGenes, rcvGenes - selected intracellular genes per side
%     sndScore, rcvScore - their selection scores (NaN if explicit)
%     sndVars, rcvVars   - [ligands, sndGenes] / [receptors, rcvGenes]
%     sndCols, rcvCols   - sndVars/rcvVars with shared names disambiguated
%     metaAll        - struct('name', 'values') for every requested
%                      Metadata attribute, full sce length
%     pivS, pivR     - metaAll restricted to one side, dropping any
%                      attribute constant there (struct('name','values'))

labels = i_grouplabels(sce, o.GroupBy);
isSnd = ismember(labels, senders);
isRcv = ismember(labels, receivers);

if ~any(isSnd)
    error('sc_sce2causalccc:noSenders', ...
        'No cell matches sender label(s) %s in "%s". Available: %s', ...
        strjoin(senders, ', '), o.GroupBy, ...
        strjoin(unique(labels)', ', '));
end
if ~any(isRcv)
    error('sc_sce2causalccc:noReceivers', ...
        'No cell matches receiver label(s) %s in "%s". Available: %s', ...
        strjoin(receivers, ', '), o.GroupBy, ...
        strjoin(unique(labels)', ', '));
end
if any(isSnd & isRcv)
    error('sc_sce2causalccc:overlap', ...
        ['%d cell(s) are labelled both sender and receiver. The two ', ...
         'populations must be disjoint.'], sum(isSnd & isRcv));
end

Xs = full(double(sce.X(:, isSnd)));
Xr = full(double(sce.X(:, isRcv)));
g = string(sce.g(:));
gU = upper(g);

if ~i_lookslikecounts(Xs)
    warning('sc_sce2causalccc:notCounts', ...
        ['sce.X does not look like raw counts (non-integer values found). ', ...
         'causalCCC asks for raw counts; export before normalising if possible.']);
end

LR = i_loadlr(o.LigandReceptor);
[pairs, ligands, receptors] = i_selectpairs(LR, gU, Xs, Xr, o.MinFrac, o.MaxPairs);
if isempty(pairs)
    error('sc_sce2causalccc:noPairs', ...
        ['No ligand-receptor pair survives: no pair has its ligand in ', ...
         '>=%.0f%% of sender cells and its receptor in >=%.0f%% of ', ...
         'receiver cells. Lower MinFrac or supply LigandReceptor.'], ...
        100*o.MinFrac, 100*o.MinFrac);
end

metaAll = i_metavalues(sce, o.Metadata);
pivS = i_metaforside(metaAll, isSnd);
pivR = i_metaforside(metaAll, isRcv);

[sndGenes, sndScore] = i_pickgenes(o.SenderGenes, Xs, g, ligands, ...
    [ligands, receptors], o.GeneSelection, o.MaxGenes, ...
    'sender', {pivS.values}, o.MIBins, o.MIPermutations);
[rcvGenes, rcvScore] = i_pickgenes(o.ReceiverGenes, Xr, g, receptors, ...
    [ligands, receptors], o.GeneSelection, o.MaxGenes, ...
    'receiver', {pivR.values}, o.MIBins, o.MIPermutations);

sndVars = [ligands, sndGenes];
rcvVars = [receptors, rcvGenes];
[sndCols, rcvCols] = i_disambiguate(sndVars, rcvVars);

ord = [find(isSnd(:)); find(isRcv(:))];

A = struct();
A.isSnd = isSnd;
A.isRcv = isRcv;
A.ord = ord;
A.g = g;
A.gU = gU;
A.Xs = Xs;
A.Xr = Xr;
A.pairs = pairs;
A.ligands = ligands;
A.receptors = receptors;
A.sndGenes = sndGenes;
A.rcvGenes = rcvGenes;
A.sndScore = sndScore;
A.rcvScore = rcvScore;
A.sndVars = sndVars;
A.rcvVars = rcvVars;
A.sndCols = sndCols;
A.rcvCols = rcvCols;
A.metaAll = metaAll;
A.pivS = pivS;
A.pivR = pivR;
end


function labels = i_grouplabels(sce, groupBy)
% Per-cell grouping labels, from the named annotation.
switch lower(groupBy)
    case {"celltype", "cell_type", "c_cell_type_tx"}
        labels = string(sce.c_cell_type_tx);
    case {"cluster", "clusterid", "c_cluster_id"}
        labels = string(sce.c_cluster_id);
    case {"batch", "batchid", "c_batch_id"}
        labels = string(sce.c_batch_id);
    otherwise
        v = sce.getCellAttribute(char(groupBy));
        if isempty(v)
            avail = string(sce.list_cell_attributes(1:2:end));
            error('sc_sce2causalccc:badGroupBy', ...
                ['GroupBy "%s" is not a known grouping. Use "celltype", ', ...
                 '"cluster", "batch", or a cell attribute: %s'], ...
                groupBy, strjoin(avail, ', '));
        end
        labels = string(v);
end
labels = labels(:);
end


function tf = i_lookslikecounts(X)
% Cheap check on a sample - a full integer test on a big matrix is wasteful
% and this only drives a warning.
v = X(X ~= 0);
if isempty(v), tf = true; return, end
if numel(v) > 20000, v = v(1:round(numel(v)/20000):end); end
tf = all(abs(v - round(v)) < 1e-9);
end


function LR = i_loadlr(userLR)
% N-by-2 string array of [ligand, receptor], uppercased for matching.
if isempty(userLR)
    pw = fileparts(fileparts(mfilename('fullpath')));   % .../private -> root
    f = fullfile(pw, 'assets', 'Ligand_Receptor', 'Ligand_Receptor_more.mat');
    if ~isfile(f)
        error('sc_sce2causalccc:noLRdb', ...
            'Bundled ligand-receptor database not found at %s.', f);
    end
    S = load(f, 'ligand', 'receptor');
    LR = [string(S.ligand(:)), string(S.receptor(:))];
elseif istable(userLR)
    LR = [string(userLR{:, 1}), string(userLR{:, 2})];
else
    LR = string(userLR);
    if size(LR, 2) ~= 2
        error('sc_sce2causalccc:badLR', ...
            'LigandReceptor must be an N-by-2 array of [ligand, receptor].');
    end
end
LR = upper(LR);
LR = LR(all(strlength(LR) > 0, 2), :);
LR = unique(LR, 'rows', 'stable');
end


function [pairs, ligands, receptors] = i_selectpairs(LR, gU, Xs, Xr, minFrac, maxPairs)
% Keep pairs whose ligand is detected in the senders and receptor in the
% receivers, then rank by the product of the two mean expressions - the
% usual first-pass CCC score, and enough to fill a MaxPairs budget.
pairs = strings(0, 2); ligands = strings(1, 0); receptors = strings(1, 0);

fracS = mean(Xs > 0, 2);
fracR = mean(Xr > 0, 2);
meanS = mean(Xs, 2);
meanR = mean(Xr, 2);

[hasL, iL] = ismember(LR(:, 1), gU);
[hasR, iR] = ismember(LR(:, 2), gU);
ok = hasL & hasR;
ok(ok) = fracS(iL(ok)) >= minFrac & fracR(iR(ok)) >= minFrac;
if ~any(ok), return, end

idx = find(ok);
score = meanS(iL(idx)) .* meanR(iR(idx));
[~, ord] = sort(score, 'descend');
idx = idx(ord);
idx = idx(1:min(maxPairs, numel(idx)));

pairs = LR(idx, :);
ligands = unique(pairs(:, 1), 'stable')';
receptors = unique(pairs(:, 2), 'stable')';
end


function [sel, score] = i_pickgenes(explicit, X, g, anchors, exclude, ...
    method, maxGenes, side, metaPivots, nbins, nperm)
% Intracellular genes for one side. An explicit list wins; otherwise rank
% by the chosen criterion and take the top maxGenes.
gU = upper(g);
score = [];

if ~isempty(explicit)
    want = upper(explicit(:)');
    missing = want(~ismember(want, gU));
    if ~isempty(missing)
        warning('sc_sce2causalccc:missingGenes', ...
            '%s gene(s) not in sce.g and skipped: %s', ...
            side, strjoin(missing, ', '));
    end
    sel = want(ismember(want, gU));
    sel = setdiff(sel, upper(exclude), 'stable');
    score = nan(1, numel(sel));
    return
end

% Never offer a gene that is silent on this side, or an anchor already in.
keep = mean(X > 0, 2) > 0 & ~ismember(gU, upper(exclude));
if ~any(keep)
    sel = strings(1, 0);
    return
end
Xk = X(keep, :);
gk = gU(keep);

switch method
    case "mi"
        [~, iA] = ismember(upper(anchors(:)), gU);
        iA = iA(iA > 0);
        [ranked, ranks] = i_rankbymi(Xk, gk, X(iA, :), metaPivots, nbins, nperm);
    case "hvg"
        T = sc_hvg(Xk, gk, true, false);
        ranked = upper(string(T.genes(:)'));
        ranks = nan(1, numel(ranked));
    case "correlation"
        [ranked, ranks] = i_rankbycorr(Xk, gk, X, gU, anchors);
end

n = min(maxGenes, numel(ranked));
sel = ranked(1:n);
score = ranks(1:n);
end


function [ranked, score] = i_rankbymi(Xk, gk, anchorX, metaPivots, nbins, nperm)
% causalCCC.MIselection in miniature: score every candidate gene by its
% mutual information with each pivot, then rank by the strongest pivot.
% Taking the max rather than the sum keeps a gene that is highly
% informative about a single receptor from being diluted by the others.
pivots = {};
for i = 1:size(anchorX, 1)
    pivots{end+1} = pkg.i_binvector(anchorX(i, :)', nbins); %#ok<AGROW>
end
for i = 1:numel(metaPivots)
    [~, ~, idx] = unique(metaPivots{i});
    pivots{end+1} = uint8(idx); %#ok<AGROW>
end

if isempty(pivots)
    ranked = gk(:)';
    score = zeros(1, numel(gk));
    return
end

Gb = i_binmatrix(Xk', nbins);
n = size(Gb, 1);

% Ranking tens of thousands of genes by MI is an extreme order statistic:
% the top score is the max over genes x pivots, so it sits well above zero
% even when nothing is dependent. On a small population that noise floor
% approaches the real signal - measured on a 56-cell population, the
% permuted null still reached ~85% of the real median.
if n < 100
    warning('sc_sce2causalccc:smallPopulation', ...
        ['MI gene selection on %d cells is close to its noise floor; the ', ...
         'ranking will favour highly detected genes. Consider pooling ', ...
         'more cells, lowering MIBins, or passing an explicit gene list.'], n);
end

% A private stream, so the calibration is reproducible without disturbing
% whatever the caller has set up.
rs = RandStream('threefry', 'Seed', 20260829);

best = zeros(size(Gb, 2), 1);
for i = 1:numel(pivots)
    mi = i_mivec(Gb, nbins, pivots{i});
    if nperm > 0
        nullmi = zeros(size(mi));
        for r = 1:nperm
            nullmi = nullmi + i_mivec(Gb, nbins, pivots{i}(randperm(rs, n)));
        end
        mi = mi - nullmi / nperm;   % excess over what chance alone gives
    end
    best = max(best, mi);
end
best = max(best, 0);

[score, ord] = sort(best, 'descend');
ranked = gk(ord)';
score = score(:)';
end


function mi = i_mivec(Gb, nb, p)
% MI, in bits, between one binned pivot and every binned candidate at once.
% The joint counts come out of a single matrix product per candidate bin,
% which is what makes pivots-against-all affordable at 20k genes.
n = size(Gb, 1);
nG = size(Gb, 2);
np = double(max(p));
if np < 2
    mi = zeros(nG, 1);   % a constant pivot informs about nothing
    return
end

P = zeros(np, n);
for b = 1:np
    P(b, :) = (p(:)' == b);
end

C = zeros(np, nb, nG);
for v = 1:nb
    C(:, v, :) = P * double(Gb == v);
end

px = sum(C, 2);                     % np x 1 x nG
py = sum(C, 1);                     % 1  x nb x nG
Cn = C / n;
denom = (px .* py) / n^2;           % implicit expansion
t = Cn .* log2(max(Cn, realmin) ./ max(denom, realmin));
t(C == 0) = 0;
mi = reshape(sum(sum(t, 1), 2), nG, 1);

% Miller-Madow: the plug-in estimate is biased upward by roughly the number
% of occupied cells, and zero-inflation leaves different genes with
% different effective support - so without this the ranking partly measures
% how many distinct values a gene happens to take.
nzJoint = reshape(sum(sum(C > 0, 1), 2), nG, 1);
nzX = reshape(sum(px > 0, 1), nG, 1);
nzY = reshape(sum(py > 0, 2), nG, 1);
mi = mi - (nzJoint - nzX - nzY + 1) / (2 * n * log(2));
mi = max(mi, 0);
end


function B = i_binmatrix(Xt, nb)
% Xt is cells-by-genes.
B = ones(size(Xt), 'uint8');
for j = 1:size(Xt, 2)
    B(:, j) = pkg.i_binvector(Xt(:, j), nb);
end
end


function [ranked, score] = i_rankbycorr(Xk, gk, X, gU, anchors)
% Rank by |correlation| with the mean profile of this side's anchors.
[hasA, iA] = ismember(upper(anchors(:)), gU);
if ~any(hasA)
    ranked = gk(:)';
    score = zeros(1, numel(gk));
    return
end
a = mean(X(iA(hasA), :), 1)';
sd = std(Xk, 0, 2);
live = sd > 0;
r = zeros(size(gk));
if std(a) > 0 && any(live)
    r(live) = abs(corr(Xk(live, :)', a));
end
r(isnan(r)) = 0;
[score, ord] = sort(r, 'descend');
ranked = gk(ord)';
score = score(:)';
end


function [sndCols, rcvCols] = i_disambiguate(sndVars, rcvVars)
% A symbol chosen for both sides needs two columns, so suffix both. Anchors
% are already excluded from the other side's pool, so this normally only
% fires for a gene that is genuinely of interest in both populations.
both = intersect(upper(sndVars), upper(rcvVars));
sndCols = sndVars;
rcvCols = rcvVars;
if isempty(both), return, end
sndCols(ismember(upper(sndVars), both)) = ...
    sndVars(ismember(upper(sndVars), both)) + "_senders";
rcvCols(ismember(upper(rcvVars), both)) = ...
    rcvVars(ismember(upper(rcvVars), both)) + "_receivers";
end


function metaAll = i_metavalues(sce, names)
% Requested metadata as full-length per-cell vectors, resolved once so both
% the MI pivots and the mosaic columns come from the same source.
metaAll = struct('name', {}, 'values', {});
for nm = names
    switch lower(nm)
        case {"batch", "batchid", "c_batch_id"}
            v = string(sce.c_batch_id);
        case {"celltype", "cell_type", "c_cell_type_tx"}
            v = string(sce.c_cell_type_tx);
        case {"cluster", "clusterid", "c_cluster_id"}
            v = string(sce.c_cluster_id);
        otherwise
            v = sce.getCellAttribute(char(nm));
            if isempty(v)
                warning('sc_sce2causalccc:missingMetadata', ...
                    'Metadata "%s" is not a cell attribute and was skipped.', nm);
                continue
            end
            v = string(v);
    end
    metaAll(end+1).name = matlab.lang.makeValidName(nm); %#ok<AGROW>
    metaAll(end).values = v(:);
end
end


function piv = i_metaforside(metaAll, mask)
% Metadata restricted to one population, dropping any that is constant
% there - a pivot with one level has zero MI with everything.
piv = struct('name', {}, 'values', {});
for k = 1:numel(metaAll)
    v = metaAll(k).values(mask);
    if numel(unique(v)) > 1
        piv(end+1).name = metaAll(k).name; %#ok<AGROW>
        piv(end).values = v;
    end
end
end


