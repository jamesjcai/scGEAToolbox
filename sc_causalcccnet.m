function [G, info] = sc_causalcccnet(sce, senders, receivers, options)
%SC_CAUSALCCCNET  Local approximation of a causalCCC/MIIC network.
%
%   G = SC_CAUSALCCCNET(sce, senders, receivers)
%   [G, info] = SC_CAUSALCCCNET(..., Name=Value)
%
%   Builds, without contacting the causalCCC web server
%   (https://miic.curie.fr/causalCCC.php) and without writing any file, a
%   local stand-in for the network that server renders from the files
%   PKG.SC_SCE2CAUSALCCC prepares. Both functions select the same senders,
%   receivers, ligand-receptor pairs and intracellular genes (via the
%   shared internal PKG.I_CAUSALCCCASSEMBLE), so this is a same-inputs sibling
%   of that exporter rather than a reimplementation of MIIC's causal
%   discovery, which this does not attempt.
%
%   What MIIC does that this does not: causal edge orientation by
%   conditional-independence testing over the whole mosaic at once,
%   including its missing-data pattern. What this does instead, and why it
%   is a reasonable stand-in:
%
%     - senders and receivers are disjoint cell populations, so no per-cell
%       statistic can relate a sender variable to a receiver variable
%       directly - the mosaic's NA pattern exists for exactly this reason.
%       Two separate intracellular networks are therefore built: one among
%       [ligands, sender genes, sender-side metadata] scored on sender
%       cells, one among [receptors, receiver genes, receiver-side
%       metadata] scored on receiver cells.
%     - the ligand-receptor pairs bridge the two: on the server they are
%       supplied as known edges, not inferred ones, and they are treated
%       the same way here - drawn directly from PKG.SC_SCE2CAUSALCCC's pair
%       selection rather than estimated from data that cannot support them.
%
%   INPUTS
%     sce, senders, receivers - as PKG.SC_SCE2CAUSALCCC
%
%   NAME-VALUE ARGUMENTS
%     GroupBy, LigandReceptor, SenderGenes, ReceiverGenes, GeneSelection,
%     MaxGenes, MaxPairs, MinFrac, Metadata, MIBins, MIPermutations
%                   as PKG.SC_SCE2CAUSALCCC; the two functions select the same
%                   variables when given the same values.
%     Method        "mi" (default) - mutual information, Miller-Madow
%                   corrected, the same estimator PKG.SC_SCE2CAUSALCCC uses for
%                   gene selection. Handles the categorical Metadata
%                   pivots directly.
%                   "pcorr"/"pearson"/"distcorr"/"xicor" - net.pcrnet,
%                   net.pearsonnet, net.distcorrnet, net.xicornet on the
%                   side's continuous expression only; Metadata is omitted
%                   for these (a warning is issued) since a categorical
%                   label has no well-defined correlation with expression.
%     TopKPerNode   keep an intracellular edge only if it is among either
%                   endpoint's TopKPerNode strongest. Default 8. An MI or
%                   correlation network among all of a side's variables is
%                   otherwise close to complete and unreadable.
%     EdgeCutoff    additionally drop intracellular edges below this
%                   |weight|, applied before TopKPerNode. Default 0
%                   (no cutoff; TopKPerNode alone does the sparsifying).
%
%   OUTPUT
%     G    - digraph. Nodes.Group is one of "ligand", "sender_gene",
%            "receptor", "receiver_gene", "metadata"; Nodes.Side is
%            "senders" or "receivers"; Nodes.Color is the hex colour
%            PKG.SC_SCE2CAUSALCCC's state-order file uses for that group, so a
%            plot coloured by Nodes.Color reads like a relative of the
%            server's own rendering. Edges.Type is "bridge" (a ligand-
%            receptor pair, direction ligand->receptor) or "inferred" (an
%            intracellular edge, direction is arbitrary - these are
%            symmetric statistics); Edges.Sign is +1/-1 for a signed
%            Method and 0 for an unsigned one (mi/distcorr/xicor).
%     info - struct with the fields PKG.SC_SCE2CAUSALCCC's info has for the
%            senders/receivers/pairs/genes it selected (nSenderCells,
%            nReceiverCells, ligands, receptors, pairs, senderGenes,
%            receiverGenes, senderGeneScores, receiverGeneScores,
%            geneSelection, groupBy), plus method (the Method used).
%
%   EXAMPLE
%     load('example_data/new_example_sce.mat', 'sce');
%     G = sc_causalcccnet(sce, "Macrophages", "Beta cells");
%     plot(G)   % a plain force-directed view; Nodes.Color/Group/Side carry
%               % the causalCCC-style grouping for a nicer rendering
%
%   For an interactive, group-coloured, sender/receiver-laid-out view, use
%   GUI.E_SCE2CAUSALCCCVIEW instead of plotting G directly.
%
%   See also: PKG.SC_SCE2CAUSALCCC, GUI.E_SCE2CAUSALCCCVIEW, NET.MINET,
%   NET.PCRNET, NET.PEARSONNET, NET.DISTCORRNET, NET.XICORNET

arguments
    sce (1,1) SingleCellExperiment
    senders (1,:) string {mustBeNonempty}
    receivers (1,:) string {mustBeNonempty}
    options.GroupBy (1,1) string = "celltype"
    options.LigandReceptor = []
    options.SenderGenes (1,:) string = strings(1, 0)
    options.ReceiverGenes (1,:) string = strings(1, 0)
    options.GeneSelection (1,1) string ...
        {mustBeMember(options.GeneSelection, ["mi", "hvg", "correlation"])} = "mi"
    options.MaxGenes (1,1) double {mustBeInteger, mustBePositive} = 15
    options.MaxPairs (1,1) double {mustBeInteger, mustBePositive} = 15
    options.MinFrac (1,1) double {mustBeNonnegative, mustBeLessThanOrEqual(options.MinFrac, 1)} = 0.10
    options.Metadata (1,:) string = strings(1, 0)
    options.MIBins (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(options.MIBins, 2)} = 5
    options.MIPermutations (1,1) double {mustBeInteger, mustBeNonnegative} = 3
    options.Method (1,1) string ...
        {mustBeMember(options.Method, ["mi", "pcorr", "pearson", "distcorr", "xicor"])} = "mi"
    options.TopKPerNode (1,1) double {mustBeInteger, mustBePositive} = 8
    options.EdgeCutoff (1,1) double {mustBeNonnegative} = 0
end

asmOpts = struct('GroupBy', options.GroupBy, 'LigandReceptor', options.LigandReceptor, ...
    'SenderGenes', options.SenderGenes, 'ReceiverGenes', options.ReceiverGenes, ...
    'GeneSelection', options.GeneSelection, 'MaxGenes', options.MaxGenes, ...
    'MaxPairs', options.MaxPairs, 'MinFrac', options.MinFrac, ...
    'Metadata', options.Metadata, 'MIBins', options.MIBins, ...
    'MIPermutations', options.MIPermutations);
A = pkg.i_causalcccassemble(sce, senders, receivers, asmOpts);

[sndNodes, sndEdges] = i_sidenetwork(A.sndVars, A.sndCols, numel(A.ligands), ...
    "ligand", "sender_gene", A.gU, A.Xs, A.pivS, "senders", ...
    options.Method, options.MIBins, options.TopKPerNode, options.EdgeCutoff);
[rcvNodes, rcvEdges] = i_sidenetwork(A.rcvVars, A.rcvCols, numel(A.receptors), ...
    "receptor", "receiver_gene", A.gU, A.Xr, A.pivR, "receivers", ...
    options.Method, options.MIBins, options.TopKPerNode, options.EdgeCutoff);

% The ligand-receptor pairs are known edges, not inferred ones - see the
% header comment - so they are added as-is rather than scored.
linkL = pkg.i_mapnames(A.pairs(:, 1), A.sndVars, A.sndCols);
linkR = pkg.i_mapnames(A.pairs(:, 2), A.rcvVars, A.rcvCols);
nBridge = numel(linkL);
bridgeEdges = table([linkL, linkR], ones(nBridge, 1), zeros(nBridge, 1), ...
    repmat("bridge", nBridge, 1), 'VariableNames', {'EndNodes', 'Weight', 'Sign', 'Type'});

nodeTable = [sndNodes; rcvNodes];
nodeTable.Color = arrayfun(@i_groupcolor, nodeTable.Group);
edgeTable = [sndEdges; rcvEdges; bridgeEdges];

G = digraph(edgeTable, nodeTable);

info = struct();
info.nSenderCells = sum(A.isSnd);
info.nReceiverCells = sum(A.isRcv);
info.ligands = A.ligands;
info.receptors = A.receptors;
info.pairs = A.pairs;
info.senderGenes = A.sndGenes;
info.receiverGenes = A.rcvGenes;
info.senderGeneScores = A.sndScore;
info.receiverGeneScores = A.rcvScore;
info.geneSelection = options.GeneSelection;
info.groupBy = options.GroupBy;
info.method = options.Method;
end


function [nodeT, edgeT] = i_sidenetwork(vars, cols, nAnchor, anchorGroup, ...
    geneGroup, gU, Xside, piv, sideTag, method, mibins, topk, cutoff)
% One side's node table (ligands|receptors + intracellular genes + this
% side's informative metadata) and its intracellular edge table.
nGene = numel(cols);
geneNodeGroup = repmat(geneGroup, nGene, 1);
geneNodeGroup(1:nAnchor) = anchorGroup;

Xg = i_extractgenerows(vars, gU, Xside);
nMeta = numel(piv);

if method == "mi"
    codes = zeros(nGene + nMeta, size(Xside, 2));
    for j = 1:nGene
        codes(j, :) = double(pkg.i_binvector(Xg(j, :)', mibins))';
    end
    metaNames = strings(nMeta, 1);
    for k = 1:nMeta
        [~, ~, idx] = unique(piv(k).values);
        codes(nGene + k, :) = idx(:)';
        metaNames(k) = string(piv(k).name) + "_" + sideTag;
    end
    names = [cols(:); metaNames];
    group = [geneNodeGroup(:); repmat("metadata", nMeta, 1)];
    W = i_pairwisemi(codes);
else
    if nMeta > 0
        warning('sc_causalcccnet:noMetadataForMethod', ...
            ['Method "%s" does not support the categorical Metadata ', ...
             'pivot(s); they are omitted from the %s-side network.'], ...
            method, sideTag);
    end
    names = cols(:);
    group = geneNodeGroup(:);
end

% Only pcorr/pearson can produce a negative weight; mi/distcorr/xicor are
% nonnegative dependence measures, so their edges carry no sign.
switch method
    case "mi"
        signed = false;
    case "pcorr"
        W = net.pcrnet(Xg);
        signed = true;
    case "pearson"
        W = full(net.pearsonnet(Xg, 0));
        signed = true;
    case "distcorr"
        W = net.distcorrnet(Xg);
        signed = false;
    case "xicor"
        W = net.xicornet(Xg, true);
        signed = false;
end

n = numel(names);
if n < 1
    nodeT = table(strings(0, 1), strings(0, 1), strings(0, 1), ...
        'VariableNames', {'Name', 'Group', 'Side'});
    edgeT = table(strings(0, 2), zeros(0, 1), zeros(0, 1), strings(0, 1), ...
        'VariableNames', {'EndNodes', 'Weight', 'Sign', 'Type'});
    return
end

W(1:n+1:end) = 0;   % no self-loops
W = i_sparsify(W, topk, cutoff);

nodeT = table(names, group, repmat(string(sideTag), n, 1), ...
    'VariableNames', {'Name', 'Group', 'Side'});
edgeT = i_edgesfromadj(W, names);
if ~signed
    edgeT.Sign(:) = 0;
end
edgeT.Type = repmat("inferred", height(edgeT), 1);
end


function W = i_sparsify(W, topk, cutoff)
n = size(W, 1);
if cutoff > 0
    W(abs(W) < cutoff) = 0;
end
if topk > 0 && topk < n - 1
    keep = false(n);
    for i = 1:n
        [~, ord] = sort(abs(W(i, :)), 'descend');
        ord = ord(1:min(topk, n-1));
        keep(i, ord(abs(W(i, ord)) > 0)) = true;
    end
    keep = keep | keep';
    W(~keep) = 0;
end
end


function T = i_edgesfromadj(W, names)
[i, j] = find(triu(W, 1));
w = zeros(numel(i), 1);
s = zeros(numel(i), 1);
for k = 1:numel(i)
    w(k) = abs(W(i(k), j(k)));
    s(k) = sign(W(i(k), j(k)));
end
EndNodes = [names(i), names(j)];
T = table(EndNodes, w, s, 'VariableNames', {'EndNodes', 'Weight', 'Sign'});
end


function c = i_groupcolor(grp)
% The exact palette PKG.SC_SCE2CAUSALCCC's state-order file assigns per group.
switch grp
    case "ligand",        c = "#2e6e5e";
    case "sender_gene",   c = "#a9dcc9";
    case "receptor",      c = "#a05a2c";
    case "receiver_gene", c = "#f3e6d8";
    case "metadata",      c = "#8b938c";
    otherwise,            c = "#999999";
end
end
