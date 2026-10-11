function [T, S] = xctnullcontrol(sce, celltype1, celltype2, opts)
%XCTNULLCONTROL  Does a scTenifoldXct result come from the gene networks?
%   T = TEN.XCTNULLCONTROL(SCE, CELLTYPE1, CELLTYPE2) runs
%   TEN.SCTENIFOLDXCT for CELLTYPE1 -> CELLTYPE2, then reruns it NumShuffles
%   times with both gene networks' labels shuffled (each side
%   independently), and scores every candidate pair's real distance against
%   its distances under the shuffles.
%
%   [T, S] = TEN.XCTNULLCONTROL(...) also returns a diagnostics struct with
%   the verdict.
%
%   WHY THIS IS NEEDED. scTenifoldXct's p-values compare candidate pairs
%   with every other gene pair. Under the default w12mode "lr" the
%   candidates are the only pairs joined by correspondence edges, so they
%   test significant whether or not the networks say anything; and the +1
%   GRN offset leaves the network signal at a few percent of the alignment
%   graph. Measured on 17 cell-type pairs in 15 datasets (pipelines/
%   glyco_xct/Notes.md section 6k): with the networks' gene labels shuffled,
%   every call was unchanged and the distances still correlated 0.997
%   (median) with the real run. A shuffled network keeps every weight and
%   every degree and only puts them on the wrong genes, so what survives a
%   shuffle is not coming from the networks.
%
%   So the number worth reporting is whether the pair beats its own
%   shuffles, and the verdict says whether the run as a whole does.
%
%   EMPP IS CALIBRATED HERE (unlike TEN.XCT2NULLCONTROL's). Checked on
%   example_data/new_example_sce.mat, Macrophages -> Endothelial cells,
%   w12mode "lr", 20 shuffles, Seed 1: with the gene NAMES permuted and the
%   expression left alone - the real networks, but the L-R list matching
%   random gene pairs with the same pair structure - five decoy sets put 9,
%   13, 15, 12 and 8 of 247 pairs at EmpP <= 0.05 (mean 11.4 against the
%   11.8 of NumExpectedByChance), pooled EmpP was flat over its 21 values
%   (chi-square p = 0.37) and the median ShuffleRatio was 1.00. The real
%   labels put 41 pairs there. The same decoys also had all 247 pairs
%   called at p_value < 0.05, as the real pairs do: under "lr" any pair
%   list joined by correspondence edges tests significant, which is why
%   p_value cannot be read as evidence and EmpP can.
%
%   INPUTS:
%     sce        - SingleCellExperiment with both cell types
%     celltype1  - source cell type
%     celltype2  - target cell type
%
%   NAME-VALUE ARGUMENTS:
%     NumShuffles - network shuffles (default 20). The empirical p cannot
%                   fall below 1/(NumShuffles+1): 20 is the least that can
%                   reach 0.05.
%     XctArgs     - cell of further name-value arguments for
%                   TEN.SCTENIFOLDXCT, used for the real run and every
%                   shuffle alike, e.g. {'w12mode', 'outer'}. 'grn1',
%                   'grn2', 'candidates' and 'verbose' are set here.
%     Grns        - the grns struct a previous TEN.SCTENIFOLDXCT call
%                   returned for the same SCE and cell types (second
%                   output). Its networks are reused rather than rebuilt.
%     Seed        - seed for the shuffles (default [], meaning the caller's
%                   random stream governs). Given, it is used through a
%                   stream of its own, so the global stream is untouched.
%     Verbose     - print progress and the verdict (default true)
%
%   OUTPUTS:
%     T - the real run's table (ligand, receptor, dist, correspondence,
%         p_value) with, per pair:
%           ShuffleMedianDist  median distance over the shuffles
%           ShuffleRatio       dist / ShuffleMedianDist: how much closer the
%                              real networks place the pair (below 1) or
%                              further (above 1). Read it beside EmpP: a
%                              small but consistent shift reaches the
%                              smallest EmpP too - in outer mode on the
%                              pancreas example, 60 pairs did while sitting
%                              only 0.2-2.8% closer than their shuffles
%           EmpP               (1 + #shuffles placing the pair at least as
%                              close) / (1 + #shuffles), one-sided
%           ShuffleSigFrac     share of shuffles in which the pair is still
%                              called at p_value < 0.05. Near 1 for a pair
%                              called in the real run, the call does not
%                              need the networks.
%         sorted by EmpP, then dist.
%     S - diagnostics struct: Verdict, NetworkSensitive, NumShuffles,
%         NumCandidates, MedianRho (Spearman of the candidates' distances,
%         real vs each shuffle, median over shuffles), MedianSigJaccard
%         (overlap of the significant sets, real vs each shuffle, median),
%         NumSigReal, NumEmpPBelow05, NumExpectedByChance (how many pairs
%         would reach EmpP <= 0.05 if the networks did nothing), MinEmpP,
%         Rho, SigJaccard.
%
%   VERDICT. NETWORK-INSENSITIVE when the shuffles leave both the ranking
%   (MedianRho >= 0.9) and the calls (MedianSigJaccard >= 0.8) in place:
%   the result then does not depend on the networks and should not be
%   reported as network-derived. Otherwise NETWORK-SENSITIVE, and EmpP says
%   which pairs it is.
%
%   One direction per call; for both, call again with the cell types
%   swapped. candidates="all" is not supported.
%
%   Example:
%       [T, grns] = ten.sctenifoldxct(sce, "Alpha cells", "Beta cells", false);
%       [Tn, S] = ten.xctnullcontrol(sce, "Alpha cells", "Beta cells", ...
%           Grns=grns, Seed=1);
%       disp(S.Verdict)
%
% see also: TEN.SCTENIFOLDXCT, TEN.KNKNULLCONTROL

arguments
    sce SingleCellExperiment
    celltype1 (1, 1) string
    celltype2 (1, 1) string
    opts.NumShuffles (1, 1) double {mustBeInteger, mustBePositive} = 20
    opts.XctArgs cell = {}
    opts.Grns struct = struct([])
    opts.Seed double {mustBeScalarOrEmpty} = []
    opts.Verbose (1, 1) logical = true
end

reserved = ["grn1", "grn2", "grn1_processed", "grn2_processed", "candidates", "verbose"];
xctArgs = opts.XctArgs;
if ~isempty(xctArgs) && any(ismember(lower(string(xctArgs(1:2:end))), reserved))
    error("TEN:XCTNULLCONTROL:ReservedArg", ...
        "XctArgs cannot set %s; this function sets them. Pass the networks " + ...
        "as Grns instead.", strjoin(reserved, ", "));
end

if isempty(opts.Seed)
    stream = RandStream.getGlobalStream();
else
    stream = RandStream("twister", "Seed", opts.Seed);
end

% -- The real run --------------------------------------------------------
if isempty(opts.Grns)
    [Treal, grns] = i_xct(sce, celltype1, celltype2, xctArgs);
else
    grns = opts.Grns;
    Treal = i_xct(sce, celltype1, celltype2, [xctArgs, {'grn1', grns.A_s, ...
        'grn1_processed', true, 'grn2', grns.A_t, 'grn2_processed', true}]);
end
nCand = height(Treal);
if nCand < 3
    error("TEN:XCTNULLCONTROL:TooFewCandidates", ...
        "Only %d candidate pair(s); a ranking comparison needs at least 3.", nCand);
end
key = Treal.ligand + "|" + Treal.receptor;
sigReal = Treal.p_value < 0.05;
ng = size(grns.A_s, 1);

% -- The shuffles --------------------------------------------------------
N = opts.NumShuffles;
D = nan(nCand, N);
sig = false(nCand, N);
rho = nan(N, 1);
jac = nan(N, 1);
if opts.Verbose
    fprintf('[xctnull] %s -> %s: %d candidate pairs, %d network shuffles ', ...
        celltype1, celltype2, nCand, N);
end
for b = 1:N
    p1 = randperm(stream, ng);
    p2 = randperm(stream, ng);
    Tb = i_xct(sce, celltype1, celltype2, [xctArgs, {'grn1', grns.A_s(p1, p1), ...
        'grn1_processed', true, 'grn2', grns.A_t(p2, p2), 'grn2_processed', true}]);
    [found, j] = ismember(key, Tb.ligand + "|" + Tb.receptor);
    D(found, b) = Tb.dist(j(found));
    sig(found, b) = Tb.p_value(j(found)) < 0.05;
    rho(b) = corr(Treal.dist(found), D(found, b), "Type", "Spearman");
    jac(b) = i_jaccard(sigReal(found), sig(found, b));
    if opts.Verbose && mod(b, 5) == 0
        fprintf('.');
    end
end
if opts.Verbose
    fprintf('\n');
end

% -- Score ---------------------------------------------------------------
asClose = sum(D <= Treal.dist, 2, "omitnan");
nScored = sum(~isnan(D), 2);
T = Treal;
T.ShuffleMedianDist = median(D, 2, "omitnan");
T.ShuffleRatio = T.dist./T.ShuffleMedianDist;
T.EmpP = (1 + asClose)./(1 + nScored);
T.ShuffleSigFrac = sum(sig, 2)./max(nScored, 1);
T = sortrows(T, ["EmpP", "dist"]);

medianRho = median(rho, "omitnan");
medianJac = median(jac, "omitnan");
insensitive = medianRho >= 0.9 && medianJac >= 0.8;
nEmp = nnz(T.EmpP <= 0.05);
% How many would get there by chance. EmpP takes the values k/(N+1), and
% under the null each is equally likely, so P(EmpP <= 0.05) is
% floor(0.05*(N+1))/(N+1): 1/21 at N = 20. Over hundreds of candidates that
% is a dozen or more pairs, so the count means something only beside it.
nChance = nCand*floor(0.05*(N + 1))/(N + 1);
chance = sprintf(' (about %.0f expected by chance)', nChance);
% Two empty call sets "agree" completely; say so rather than quote a Jaccard
if ~any(sigReal) && ~any(sig, "all")
    calls = 'no pair is called at p < 0.05 in the real run or any shuffle';
else
    calls = sprintf('the significant calls overlap %.0f%% (median Jaccard)', ...
        100*medianJac);
end
if insensitive
    verdict = sprintf(['NETWORK-INSENSITIVE - with the networks'' gene labels ' ...
        'shuffled, the distances still correlate %.3f (median Spearman) with ' ...
        'the real run, and %s. This result does not depend on the networks; ' ...
        'do not report it as network-derived. %d of %d pairs beat their ' ...
        'shuffles at EmpP <= 0.05%s.'], medianRho, calls, nEmp, nCand, chance);
else
    verdict = sprintf(['NETWORK-SENSITIVE - shuffling the networks moves the ' ...
        'distances (median Spearman %.3f with the real run); %s. %d of %d ' ...
        'pairs beat their shuffles at EmpP <= 0.05%s; read EmpP pair by pair ' ...
        'only beyond that.'], medianRho, calls, nEmp, nCand, chance);
end
if N < 19
    verdict = [verdict, sprintf([' (With %d shuffles EmpP cannot go below ' ...
        '%.3f; use NumShuffles >= 19 to read it at 0.05.)'], N, 1/(N + 1))];
end
if opts.Verbose
    fprintf('[xctnull] %s\n', verdict);
end

S = struct("Verdict", string(verdict), "NetworkSensitive", ~insensitive, ...
    "NumShuffles", N, "NumCandidates", nCand, "MedianRho", medianRho, ...
    "MedianSigJaccard", medianJac, "NumSigReal", nnz(sigReal), ...
    "NumEmpPBelow05", nEmp, "NumExpectedByChance", nChance, ...
    "MinEmpP", 1/(N + 1), "Rho", rho, "SigJaccard", jac);
end


function [T, grns] = i_xct(sce, a, b, args) %#ok<INUSD> used inside evalc
% One direction, database candidates, quietly
[~, T, grns] = evalc("ten.sctenifoldxct(sce, a, b, false, 'verbose', false, args{:})");
end


function j = i_jaccard(x, y)
% Overlap of two logical sets; two empty sets agree completely
u = nnz(x | y);
if u == 0
    j = 1;
else
    j = nnz(x & y)/u;
end
end
