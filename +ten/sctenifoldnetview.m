function [hFig, info] = sctenifoldnetview(net1, net2, genes, T, options)
% SCTENIFOLDNETVIEW  Side-by-side viewer for a scTenifoldNet result.
%
%   ten.sctenifoldnetview(net1, net2) draws the two GRNs side by side over a
%   shared node layout, so a node sits in the same place in both panels and
%   the two can be read position by position. Set ShowDifference to add a
%   third panel holding net2 minus net1. NET1 and NET2 are either NxN
%   adjacency matrices or paths to .mat files holding a variable A (the
%   adjacency matrix) and genes (the Nx1 gene list) - the format written by
%   ten.sctenifoldnet and gui.callback_scTenifoldNet1lite.
%
%   ten.sctenifoldnetview(net1, net2, genes) supplies the gene list
%   explicitly. Required when net1/net2 are matrices, optional when they
%   are .mat paths.
%
%   ten.sctenifoldnetview(net1, net2, genes, T) uses the differential
%   regulation table T from ten.sctenifoldnet / ten.i_dr - columns
%   genelist, drdist, FC, pValues, pAdjusted - to pick and score the seed
%   genes. T may also be a path to a .csv holding that table. With T
%   omitted or empty, genes are ranked by how much their edges moved
%   between the two networks instead.
%
%   3000 genes cannot be drawn legibly, so the viewer shows a neighbourhood
%   rather than the whole network: the top differentially regulated genes
%   are taken as seeds, then their strongest partners in either network are
%   added until NumNodes nodes are reached.
%
%   Name-value arguments
%     Mode         "neighborhood" (default) seeds with the top NumSeeds DR
%                  genes and grows outwards. "drmodule" instead takes EVERY
%                  gene passing PAdjCutoff, so the figure is the whole
%                  differentially regulated set in context rather than a
%                  top-N slice of it; NumSeeds is then ignored, and NumNodes
%                  only governs how many partners are added on top - a
%                  significant gene is never dropped to meet it.
%     SeedGenes    Seed genes, overriding the top-DR selection.
%     NumSeeds     Number of top-DR genes used as seeds (default 10). Only
%                  genes with pAdjusted below PAdjCutoff are eligible; if
%                  fewer qualify, the top NumSeeds by drdist are used.
%                  Ignored when Mode is "drmodule".
%     NumNodes     Target node count after expansion (default 60).
%     PAdjCutoff   Adjusted p-value defining a significant DR gene (0.05).
%     ExpandBy     How partners are ranked when the seeds are expanded:
%                  "strength" (default) takes each seed's strongest
%                  partners in either condition, which shows the seed in
%                  its regulatory context; "change" takes the partners
%                  whose edge to the seed moved the most.
%     EdgeCutoff   Quantile in [0,1) of |weight| below which edges are
%                  hidden (default 0.9, i.e. draw the strongest 10%).
%                  A single threshold is shared by the two condition
%                  panels so their edge counts are comparable.
%     MinDegree    Every node keeps its MinDegree strongest edges even if
%                  they fall under EdgeCutoff (default 2). Without this a
%                  DR gene is routinely drawn with no edges at all: it was
%                  selected for how much its edges moved, which says
%                  nothing about how strong they are.
%     ShowDifference  Add a third panel holding net2 minus net1 (default
%                  false). It shares the node set and layout of the other
%                  two but not their edge threshold: the differences are on
%                  a smaller scale and get their own quantile. A node with
%                  no edges there is a gene whose edges did not move.
%     Labels       1x2 string naming the two conditions. Defaults to the
%                  genotype field of the .mat files when available.
%     Title        Figure title. Defaults to the celltype field.
%     Symmetrize   Replace A with (A + A')/2 before anything else (default
%                  true). Bootstrapped scTenifoldNet output is near- but not
%                  exactly symmetric, and an undirected graph is what the DR
%                  test assumes. The panels are drawn undirected either way;
%                  turning this off only changes which weight is plotted.
%     ParentFig    Parent figure for theme and placement.
%
%   [hFig, info] = ... also returns a struct with the node set, the seed
%   genes, the subnetwork adjacency matrices, the graph objects and the
%   thresholds actually used. info.Graphs and info.EdgeThreshold have one
%   entry per panel, so two unless ShowDifference is set; info.Adiff is
%   empty in that case. info.WeaklyConnectedSeeds names the seeds that
%   ended up with no edge in either network at the chosen EdgeCutoff.
%
%   Example
%     ten.sctenifoldnetview('grn_Microglia_WT.mat', 'grn_Microglia_KO.mat', ...
%         [], 'dr_Microglia.csv', 'NumSeeds', 8, 'NumNodes', 60);
%
%     % with the difference panel
%     ten.sctenifoldnetview('grn_Microglia_WT.mat', 'grn_Microglia_KO.mat', ...
%         [], 'dr_Microglia.csv', 'ShowDifference', true);
%
%     % every significant DR gene, and nothing else
%     ten.sctenifoldnetview('grn_Microglia_WT.mat', 'grn_Microglia_KO.mat', ...
%         [], 'dr_Microglia.csv', 'Mode', 'drmodule', 'NumNodes', 1);
%
% See also ten.sctenifoldnet, ten.i_dr, sc_grnview2, gui.i_multigraphs.

arguments
    net1
    net2
    genes = []
    T = []
    options.SeedGenes string = string.empty
    options.NumSeeds (1, 1) double {mustBePositive, mustBeInteger} = 10
    options.NumNodes (1, 1) double {mustBePositive, mustBeInteger} = 60
    options.PAdjCutoff (1, 1) double {mustBePositive} = 0.05
    options.EdgeCutoff (1, 1) double {mustBeGreaterThanOrEqual(options.EdgeCutoff, 0), mustBeLessThan(options.EdgeCutoff, 1)} = 0.9
    options.MinDegree (1, 1) double {mustBeNonnegative, mustBeInteger} = 2
    options.ShowDifference (1, 1) logical = false
    options.Mode (1, 1) string {mustBeMember(options.Mode, ["neighborhood", "drmodule"])} = "neighborhood"
    options.ExpandBy (1, 1) string {mustBeMember(options.ExpandBy, ["strength", "change"])} = "strength"
    options.Labels string = string.empty
    options.Title string = string.empty
    options.Symmetrize (1, 1) logical = true
    options.ParentFig = []
end

[A0, g0, meta0] = i_resolvenet(net1, 'net1');
[A1, g1, meta1] = i_resolvenet(net2, 'net2');

if isempty(genes)
    genes = g0;
    if isempty(genes), genes = g1; end
end
if isempty(genes)
    error('MATLAB:sctenifoldnetview:noGeneList', ...
        'No gene list found. Pass genes, or use .mat files that store one.');
end
genes = string(genes(:));

n = size(A0, 1);
if size(A0, 2) ~= n || ~isequal(size(A1), size(A0))
    error('MATLAB:sctenifoldnetview:sizeMismatch', ...
        'The two adjacency matrices must be square and the same size.');
end
if numel(genes) ~= n
    error('MATLAB:sctenifoldnetview:geneCount', ...
        'genes has %d entries but the networks are %dx%d.', numel(genes), n, n);
end
if ~isempty(g0) && ~isempty(g1) && ~isequal(g0, g1)
    error('MATLAB:sctenifoldnetview:geneListMismatch', ...
        'The two networks are indexed by different gene lists.');
end

A0 = double(A0);
A1 = double(A1);
if options.Symmetrize
    A0 = 0.5 * (A0 + A0.');
    A1 = 0.5 * (A1 + A1.');
end

labels = options.Labels;
if isempty(labels)
    labels = [i_metastr(meta0, 'genotype', "Network 1"), ...
        i_metastr(meta1, 'genotype', "Network 2")];
end
labels = string(labels(:))';
if numel(labels) ~= 2
    error('MATLAB:sctenifoldnetview:badLabels', ...
        'Labels must name exactly two conditions.');
end

figtitle = options.Title;
if isempty(figtitle)
    figtitle = i_metastr(meta0, 'celltype', "");
end
figtitle = char(strjoin(string(figtitle), ' '));

% --- rank genes by differential regulation ------------------------------
[T, score, scoreshort] = i_resolvedrtable(T, genes, A0, A1);

% --- choose seeds -------------------------------------------------------
if isempty(options.SeedGenes) && options.Mode == "drmodule"
    seedidx = i_picksigseeds(T, genes, options.PAdjCutoff, options.NumSeeds);
elseif ~isempty(options.SeedGenes)
    seedgenes = string(options.SeedGenes(:));
    missing = ~ismember(seedgenes, genes);
    if all(missing)
        error('MATLAB:sctenifoldnetview:noValidSeeds', ...
            'None of the SeedGenes are in the gene list.');
    end
    if any(missing)
        warning('scGEAToolbox:sctenifoldnetview:missingSeeds', ...
            'Seed gene(s) not in the network, dropped: %s', ...
            strjoin(seedgenes(missing), ', '));
        seedgenes = seedgenes(~missing);
    end
    [~, seedidx] = ismember(seedgenes, genes);
else
    seedidx = i_pickseeds(T, genes, score, options.NumSeeds, options.PAdjCutoff);
end
seedidx = unique(seedidx(:), 'stable');

% --- expand to a drawable neighbourhood ---------------------------------
% In drmodule mode every significant gene stays in, so the node budget can
% only ever grow to fit them: NumNodes then governs how many partners get
% added on top, not whether a DR gene is dropped.
targetnodes = max(options.NumNodes, numel(seedidx));
[nodeidx, edgefloor] = i_expand(A0, A1, seedidx, targetnodes, options.ExpandBy);

sub0 = A0(nodeidx, nodeidx);
sub1 = A1(nodeidx, nodeidx);
subd = sub1 - sub0;
subgenes = genes(nodeidx);

% --- prune edges for display --------------------------------------------
% One threshold for both condition panels, so a panel with fewer edges
% really does have weaker edges rather than a gentler cutoff. The
% difference panel is on a different scale and gets its own.
tCond = i_pooledquantile({sub0, sub1}, options.EdgeCutoff);
sub0 = i_prune(sub0, tCond, options.MinDegree);
sub1 = i_prune(sub1, tCond, options.MinDegree);

G0 = graph(sub0, cellstr(subgenes), 'omitselfloops');
G1 = graph(sub1, cellstr(subgenes), 'omitselfloops');

% A seed can still land bare when every one of its edges is zero, which
% MinDegree cannot fix - there is nothing to keep. That is a real property
% of the gene, not a display artefact: it was selected for how far its
% edges moved, and a large relative move on a weak edge is still weak. Say
% so rather than let the reader assume the layout dropped it.
isseed = ismember(nodeidx, seedidx);
weakseeds = subgenes(isseed & degree(G0) == 0 & degree(G1) == 0);
if ~isempty(weakseeds)
    warning('scGEAToolbox:sctenifoldnetview:weaklyConnectedSeeds', ...
        ['%d seed gene(s) have no edge in either network at this cutoff ', ...
        'and are drawn unconnected: %s. Lower EdgeCutoff to reach them.'], ...
        numel(weakseeds), strjoin(weakseeds, ', '));
end

condref = max(abs([sub0(:); sub1(:)]));
Gs = {G0, G1};
paneltitles = [ ...
    string(sprintf('%s  (%d edges)', labels(1), G0.numedges)), ...
    string(sprintf('%s  (%d edges)', labels(2), G1.numedges))];
edgeref = [condref, condref];
thresholds = [tCond, tCond];

if options.ShowDifference
    tDiff = i_pooledquantile({subd}, options.EdgeCutoff);
    subd = i_prune(subd, tDiff, options.MinDegree);
    Gd = graph(subd, cellstr(subgenes), 'omitselfloops');
    Gs{end+1} = Gd;
    paneltitles(end+1) = string(sprintf('%s - %s  (%d edges)', ...
        labels(2), labels(1), Gd.numedges));
    edgeref(end+1) = max(abs(subd(:)));
    thresholds(end+1) = tDiff;
else
    subd = [];
    tDiff = [];
end

% --- decorate and draw ---------------------------------------------------
nodeinfo = struct;
nodeinfo.PanelTitles = paneltitles;
% Red marks significance, not seed status. Seeds are how the node set was
% chosen - an implementation detail - whereas significance is the claim the
% figure is making. Marking seeds also missed the case that matters most:
% a partner gene pulled in for context that is itself significant.
[nodeinfo.Highlight, nodeinfo.HighlightName, nodeinfo.PlainName] = ...
    i_highlightmask(T, subgenes, isseed, options.PAdjCutoff, scoreshort);
nodeinfo.IsTF = i_istf(subgenes);
nodeinfo.EdgeRef = edgeref;
nodeinfo.Legend = i_legendtext(options.ShowDifference, labels, ...
    tCond, tDiff, options.MinDegree);

if isempty(figtitle)
    figname = sprintf('scTenifoldNet: %s vs %s', labels(1), labels(2));
else
    figname = sprintf('scTenifoldNet: %s vs %s (%s)', labels(1), labels(2), figtitle);
end

[hFig, pv] = gui.i_multigraphs(Gs, nodeinfo, figname, options.ParentFig);

if nargout > 1
    info = struct;
    info.Genes = subgenes;
    info.NodeIdx = nodeidx;
    info.SeedGenes = genes(seedidx);
    % Still returned even though nothing on the figure encodes it now: it
    % is what ranked the genes, so a caller tabulating or re-plotting the
    % result needs it.
    info.Score = score(nodeidx);
    info.A1 = sub0;
    info.A2 = sub1;
    info.Adiff = subd;
    info.Graphs = Gs;
    info.Labels = labels;
    info.EdgeThreshold = thresholds;
    info.NeighborFloor = edgefloor;
    info.WeaklyConnectedSeeds = weakseeds;
    info.DRTable = T;
    info.PlotHandles = pv;
end
end

% ------------------------------------------------------- local functions

function [A, g, meta] = i_resolvenet(net, argname)
% Accept an adjacency matrix or a path to a .mat file holding A and genes.
g = string.empty;
meta = struct;
if isnumeric(net) || islogical(net) || issparse(net)
    A = net;
    return;
end
if ~(ischar(net) || (isstring(net) && isscalar(net)))
    error('MATLAB:sctenifoldnetview:badNetArg', ...
        '%s must be a numeric matrix or a path to a .mat file.', argname);
end
fname = char(net);
if ~isfile(fname)
    error('MATLAB:sctenifoldnetview:fileNotFound', ...
        '%s: file not found: %s', argname, fname);
end
meta = load(fname);
if ~isfield(meta, 'A')
    error('MATLAB:sctenifoldnetview:noAdjacency', ...
        '%s: %s has no variable named A.', argname, fname);
end
A = meta.A;
for name = ["genes", "g", "genelist"]
    if isfield(meta, name)
        g = string(meta.(name));
        g = g(:);
        break;
    end
end
end

function s = i_metastr(meta, field, fallback)
s = fallback;
if isstruct(meta) && isfield(meta, field) && ~isempty(meta.(field))
    s = string(meta.(field));
    s = s(1);
end
end

function [T, score, scoreshort] = i_resolvedrtable(T, genes, A0, A1)
% Return the DR table aligned to nothing in particular (it keeps its own
% gene column) plus a per-gene score in the order of genes.
n = numel(genes);
if ~isempty(T) && (ischar(T) || (isstring(T) && isscalar(T)))
    fname = char(T);
    if ~isfile(fname)
        error('MATLAB:sctenifoldnetview:fileNotFound', ...
            'DR table file not found: %s', fname);
    end
    T = readtable(fname);
end
if isempty(T)
    % No DR test available: fall back to how far each gene's edge profile
    % moved. Same quantity ten.i_dr squares, before the chi-square step.
    T = table.empty;
    score = vecnorm(A1-A0, 2, 2);
    scoreshort = "edge shift";
    return;
end
if ~istable(T)
    error('MATLAB:sctenifoldnetview:badDRTable', ...
        'T must be a table or a path to a .csv file.');
end
if ~ismember('genelist', T.Properties.VariableNames)
    error('MATLAB:sctenifoldnetview:noGenelistColumn', ...
        'The DR table has no genelist column.');
end
tgenes = string(T.genelist);
score = zeros(n, 1);
[found, loc] = ismember(genes, tgenes);
if ismember('pAdjusted', T.Properties.VariableNames)
    padj = T.pAdjusted(loc(found));
    % -log10(padj) saturates at the smallest representable p; cap so a
    % single astronomically significant gene does not flatten the rest.
    v = -log10(max(padj, realmin));
    score(found) = min(v, 50);
    scoreshort = "DR significance";
elseif ismember('drdist', T.Properties.VariableNames)
    score(found) = T.drdist(loc(found));
    scoreshort = "DR distance";
else
    error('MATLAB:sctenifoldnetview:noScoreColumn', ...
        'The DR table has neither a pAdjusted nor a drdist column.');
end
end

function seedidx = i_picksigseeds(T, genes, padjcutoff, numseeds)
% Every gene that clears the significance cutoff, best first - the whole
% set, not a top-N slice of it. Ranking by drdist and filtering by
% pAdjusted agree by construction: ten.i_dr derives the p-value from
% drdist, so the two orderings are monotonically related.
if isempty(T)
    error('MATLAB:sctenifoldnetview:drmoduleNeedsTable', ...
        ['Mode "drmodule" needs a DR table to know which genes are ', ...
        'significant. Pass T, or use the default Mode.']);
end
if ~ismember('pAdjusted', T.Properties.VariableNames)
    error('MATLAB:sctenifoldnetview:noPAdjusted', ...
        'Mode "drmodule" needs a pAdjusted column in the DR table.');
end
if ismember('drdist', T.Properties.VariableNames)
    T = sortrows(T, 'drdist', 'descend');
else
    T = sortrows(T, 'pAdjusted', 'ascend');
end
tgenes = string(T.genelist);
keep = ismember(tgenes, genes) & T.pAdjusted < padjcutoff;
if ~any(keep)
    % Fall back to the top of the ranking rather than refusing to draw.
    % Taken inline, not by calling i_pickseeds, which would re-apply the
    % same cutoff and warn a second time about the same thing.
    warning('scGEAToolbox:sctenifoldnetview:noSignificantDR', ...
        'No gene reaches pAdjusted < %g; using the top %d by rank.', ...
        padjcutoff, numseeds);
    inlist = find(ismember(tgenes, genes), numseeds);
    [~, seedidx] = ismember(tgenes(inlist), genes);
    return;
end
[~, seedidx] = ismember(tgenes(keep), genes);
end

function seedidx = i_pickseeds(T, genes, score, numseeds, padjcutoff)
if isempty(T)
    % No DR test: fall back to the same score that sizes the nodes.
    [~, order] = sort(score, 'descend');
    seedidx = order(1:min(numseeds, numel(order)));
    return;
end
if ismember('drdist', T.Properties.VariableNames)
    T = sortrows(T, 'drdist', 'descend');
elseif ismember('pAdjusted', T.Properties.VariableNames)
    T = sortrows(T, 'pAdjusted', 'ascend');
end
tgenes = string(T.genelist);
keep = ismember(tgenes, genes);
if ismember('pAdjusted', T.Properties.VariableNames)
    sig = keep & T.pAdjusted < padjcutoff;
    if sum(sig) >= 1
        keep = sig;
    else
        warning('scGEAToolbox:sctenifoldnetview:noSignificantDR', ...
            'No gene reaches pAdjusted < %g; using the top %d by rank.', ...
            padjcutoff, numseeds);
    end
end
tgenes = tgenes(keep);
tgenes = tgenes(1:min(numseeds, numel(tgenes)));
if isempty(tgenes)
    error('MATLAB:sctenifoldnetview:noSeedInNetwork', ...
        'No DR gene from the table is present in the network gene list.');
end
[~, seedidx] = ismember(tgenes, genes);
end

function [nodeidx, edgefloor] = i_expand(A0, A1, seedidx, numnodes, mode)
% Take partners one rank at a time from each seed in turn, so every seed is
% represented. Ranking whole columns at once, rather than sweeping a
% threshold, keeps the node count exact and avoids a quantile over a
% 9-million-element matrix. A plain global top-N does not work here: one
% strongly connected seed takes every slot and the rest are drawn bare.
if mode == "change"
    B = abs(A1(:, seedidx)-A0(:, seedidx));
else
    B = max(abs(A0(:, seedidx)), abs(A1(:, seedidx)));
end
B(seedidx, :) = 0; % seeds are already in
nseed = numel(seedidx);
[sortedval, sortedidx] = sort(B, 1, 'descend');

nodeidx = seedidx(:);
taken = false(size(A0, 1), 1);
taken(seedidx) = true;
edgefloor = inf;
for rank = 1:size(sortedidx, 1)
    if numel(nodeidx) >= numnodes, break; end
    for s = 1:nseed
        if sortedval(rank, s) <= 0, continue; end
        cand = sortedidx(rank, s);
        if taken(cand), continue; end
        taken(cand) = true;
        nodeidx(end+1, 1) = cand; %#ok<AGROW>
        edgefloor = min(edgefloor, sortedval(rank, s));
        if numel(nodeidx) >= numnodes, break; end
    end
end
if ~isfinite(edgefloor), edgefloor = 0; end
end

function t = i_pooledquantile(Ms, q)
% Threshold on |weight|, pooled over the off-diagonal entries of every
% matrix in Ms so that panels sharing the threshold stay comparable.
a = [];
for k = 1:numel(Ms)
    M = Ms{k};
    a = [a; abs(nonzeros(triu(M, 1)))]; %#ok<AGROW>
end
if isempty(a)
    t = 0;
else
    t = quantile(a, q);
end
end

function M = i_prune(M, t, mindeg)
% Hide edges under t, but never strip a node bare: each one keeps its
% mindeg strongest edges regardless. A DR gene is chosen for how far its
% edges moved, which is unrelated to how strong they are, so a plain
% threshold leaves the very genes the plot is about as isolated dots.
M = M - diag(diag(M));
M = 0.5 * (M + M.'); % kill the round-off asymmetry graph() rejects
keep = abs(M) >= t & M ~= 0;
if mindeg > 0
    [~, order] = sort(abs(M), 2, 'descend');
    for k = 1:size(M, 1)
        top = order(k, 1:min(mindeg, size(M, 2)));
        top = top(M(k, top) ~= 0);
        keep(k, top) = true;
        keep(top, k) = true;
    end
end
M = M .* keep;
end

function txt = i_legendtext(showdiff, labels, tCond, tDiff, mindeg)
if showdiff
    txt = sprintf(['\nIn the difference panel the weight is %s minus %s, ', ...
        'so a blue edge became more positive in %s and a red edge became ', ...
        'more negative in %s. Either can mean a stronger interaction: ', ...
        'read the sign in panels 1 and 2 to tell which.\n', ...
        'Panels 1 and 2 share an edge threshold of %.3g; the difference ', ...
        'panel uses %.3g. Every node keeps its %d strongest edges even ', ...
        'below the threshold.'], labels(2), labels(1), labels(2), ...
        labels(2), tCond, tDiff, mindeg);
else
    txt = sprintf(['\nBoth panels share an edge threshold of %.3g, so a ', ...
        'panel with fewer edges really does have weaker edges. Every node ', ...
        'keeps its %d strongest edges even below the threshold.\n', ...
        'Pass ShowDifference = true for a third panel holding %s minus %s.'], ...
        tCond, mindeg, labels(2), labels(1));
end
end

function [mask, hname, pname] = i_highlightmask(T, subgenes, isseed, padjcutoff, scoreshort)
% Prefer significance. Without p-values there is nothing to be significant
% about, so fall back to marking the seeds and say so in the legend rather
% than leaving the reader to guess what the red means.
if ~isempty(T) && ismember('pAdjusted', T.Properties.VariableNames)
    tgenes = string(T.genelist);
    [found, loc] = ismember(subgenes, tgenes);
    mask = false(numel(subgenes), 1);
    mask(found) = T.pAdjusted(loc(found)) < padjcutoff;
    hname = sprintf('significant (p.adj < %g)', padjcutoff);
    pname = "not significant";
else
    mask = isseed(:);
    hname = sprintf('seed genes (top %s)', scoreshort);
    pname = "partner genes";
end
end

function tf = i_istf(genelist)
tf = false(numel(genelist), 1);
mfolder = fileparts(mfilename('fullpath'));
fname = fullfile(mfolder, '..', 'assets', 'TFome', 'tfome_tfgenes.mat');
if ~isfile(fname), return; end
S = load(fname, 'tfgenes');
tf = ismember(upper(string(genelist)), string(S.tfgenes));
end
