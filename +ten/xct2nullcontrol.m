function [T, S] = xct2nullcontrol(sce1, sce2, celltype1, celltype2, opts)
%XCT2NULLCONTROL  Does a scTenifoldXct2 result come from the gene networks?
%   T = TEN.XCT2NULLCONTROL(SCE1, SCE2, CELLTYPE1, CELLTYPE2) runs
%   TEN.SCTENIFOLDXCT2 for CELLTYPE1 -> CELLTYPE2, then reruns it
%   NumShuffles times with the networks' gene labels shuffled, and scores
%   every candidate pair's real change between samples (diff2) against its
%   changes under the shuffles.
%
%   [T, S] = TEN.XCT2NULLCONTROL(...) also returns a diagnostics struct with
%   the verdict.
%
%   The differential counterpart of TEN.XCTNULLCONTROL; see that function
%   for why the check is needed. Each shuffle uses ONE permutation for
%   CELLTYPE1's networks in both samples and one for CELLTYPE2's. The
%   difference between the samples' networks is therefore kept, edge for
%   edge, and only put on the wrong genes. Shuffling each sample's network
%   independently would instead invent differences everywhere and make any
%   real result look specific.
%
%   INPUTS:
%     sce1, sce2  - SingleCellExperiment objects for the two conditions,
%                   sharing one gene list
%     celltype1   - source cell type
%     celltype2   - target cell type
%
%   NAME-VALUE ARGUMENTS:
%     NumShuffles - network shuffles (default 20). The empirical p cannot
%                   fall below 1/(NumShuffles+1): 20 is the least that can
%                   reach 0.05.
%     XctArgs     - cell of further name-value arguments for
%                   TEN.SCTENIFOLDXCT2, used for the real run and every
%                   shuffle alike, e.g. {'w12mode', 'outer'}. The grn_*
%                   arguments and 'candidates' are set here.
%     Grns        - the grns struct a previous TEN.SCTENIFOLDXCT2 call
%                   returned for the same samples and cell types, in this
%                   direction (A_s1/A_s2 are CELLTYPE1's networks). Reused
%                   rather than rebuilt.
%     Seed        - seed for the shuffles (default [], meaning the caller's
%                   random stream governs). Given, it is used through a
%                   stream of its own, so the global stream is untouched.
%     Verbose     - print progress and the verdict (default true)
%
%   OUTPUTS:
%     T - the real run's table (ligand, receptor, dist_s1, dist_s2, diff2,
%         FC, p_value, q_value) with, per pair:
%           ShuffleMedianDiff2  median diff2 over the shuffles
%           ShuffleRatio        diff2 / ShuffleMedianDiff2: how much more
%                               (above 1) or less the real networks make
%                               the pair change between samples
%           EmpP                (1 + #shuffles changing the pair at least as
%                               much) / (1 + #shuffles), one-sided
%           ShuffleSigFrac      share of shuffles in which the pair is still
%                               called at q_value < 0.05, the threshold
%                               llm.run_sctenifoldxct2 reports
%         sorted by EmpP, then diff2 descending.
%     S - diagnostics struct: Verdict, NetworkSensitive, NumShuffles,
%         NumCandidates, MedianRho (Spearman of the candidates' diff2, real
%         vs each shuffle, median), MedianSigJaccard (overlap of the
%         q < 0.05 sets, median), NumSigReal, NumEmpPBelow05,
%         NumExpectedByChance (pairs at EmpP <= 0.05 if real and shuffled
%         networks were exchangeable, which random splits show they are
%         not; see below), MinEmpP, Rho,
%         SigJaccard.
%
%   VERDICT. NETWORK-INSENSITIVE when the shuffles leave both the ranking
%   (MedianRho >= 0.9) and the calls (MedianSigJaccard >= 0.8) in place;
%   otherwise NETWORK-SENSITIVE.
%
%   EMPP IS NOT A TEST OF A CHANGE BETWEEN CONDITIONS. Real and shuffled
%   networks are not exchangeable: the real ones put their structure on the
%   genes the ligand-receptor links join. On six random halves of one sample
%   (example_data/new_example_sce.mat, Macrophages -> Endothelial cells,
%   20 shuffles), where nothing differs, 17-43 pairs (mean 25) reached
%   EmpP <= 0.05 against the 11.8 an exchangeable null gives
%   (NumExpectedByChance); pooled EmpP fell from 127 at 1/21 to about 40 at 1
%   where a flat 59 was expected; the real diff2 ran 1.2-3.4 times its
%   shuffled median; and every split was rated NETWORK-SENSITIVE. So
%   NETWORK-SENSITIVE says the networks move the ranking, not that the
%   conditions differ, and a count of pairs at small EmpP means something
%   only against random splits of the same cells - TEN.XCT2SPLITCALIBRATION
%   provides them. For calls, use scTenifoldXct2's own q_value, which
%   called nothing in any of those splits.
%
%   One direction per call; for both, call again with the cell types
%   swapped and Grns rearranged to match (A_s* <-> A_t*).
%   candidates="all" is not supported.
%
%   Example:
%       [T, grns] = ten.sctenifoldxct2(sce1, sce2, "Alpha cells", "Beta cells", false);
%       [Tn, S] = ten.xct2nullcontrol(sce1, sce2, "Alpha cells", "Beta cells", ...
%           Grns=grns, Seed=1);
%       disp(S.Verdict)
%
% see also: TEN.SCTENIFOLDXCT2, TEN.XCTNULLCONTROL

arguments
    sce1 SingleCellExperiment
    sce2 SingleCellExperiment
    celltype1 (1, 1) string
    celltype2 (1, 1) string
    opts.NumShuffles (1, 1) double {mustBeInteger, mustBePositive} = 20
    opts.XctArgs cell = {}
    opts.Grns struct = struct([])
    opts.Seed double {mustBeScalarOrEmpty} = []
    opts.Verbose (1, 1) logical = true
end

reserved = ["grn_s1", "grn_t1", "grn_s2", "grn_t2", "grn_s1_processed", ...
    "grn_t1_processed", "grn_s2_processed", "grn_t2_processed", "candidates"];
xctArgs = opts.XctArgs;
if ~isempty(xctArgs) && any(ismember(lower(string(xctArgs(1:2:end))), reserved))
    error("TEN:XCT2NULLCONTROL:ReservedArg", ...
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
    [Treal, grns] = i_xct2(sce1, sce2, celltype1, celltype2, xctArgs);
else
    grns = opts.Grns;
    Treal = i_xct2(sce1, sce2, celltype1, celltype2, [xctArgs, ...
        i_grnargs(grns, 1:size(grns.A_s1, 1), 1:size(grns.A_s1, 1))]);
end
if isempty(grns)
    error("TEN:XCT2NULLCONTROL:NoNetworks", ...
        "TEN.SCTENIFOLDXCT2 stopped before building networks (no L-R pairs?).");
end
nCand = height(Treal);
if nCand < 3
    error("TEN:XCT2NULLCONTROL:TooFewCandidates", ...
        "Only %d candidate pair(s); a ranking comparison needs at least 3.", nCand);
end
key = Treal.ligand + "|" + Treal.receptor;
sigReal = Treal.q_value < 0.05;
ng = size(grns.A_s1, 1);

% -- The shuffles --------------------------------------------------------
N = opts.NumShuffles;
D = nan(nCand, N);
sig = false(nCand, N);
rho = nan(N, 1);
jac = nan(N, 1);
if opts.Verbose
    fprintf('[xct2null] %s -> %s: %d candidate pairs, %d network shuffles ', ...
        celltype1, celltype2, nCand, N);
end
for b = 1:N
    pS = randperm(stream, ng);
    pT = randperm(stream, ng);
    Tb = i_xct2(sce1, sce2, celltype1, celltype2, [xctArgs, i_grnargs(grns, pS, pT)]);
    [found, j] = ismember(key, Tb.ligand + "|" + Tb.receptor);
    D(found, b) = Tb.diff2(j(found));
    sig(found, b) = Tb.q_value(j(found)) < 0.05;
    rho(b) = corr(Treal.diff2(found), D(found, b), "Type", "Spearman");
    jac(b) = i_jaccard(sigReal(found), sig(found, b));
    if opts.Verbose && mod(b, 5) == 0
        fprintf('.');
    end
end
if opts.Verbose
    fprintf('\n');
end

% -- Score ---------------------------------------------------------------
atLeast = sum(D >= Treal.diff2, 2, "omitnan");
nScored = sum(~isnan(D), 2);
T = Treal;
T.ShuffleMedianDiff2 = median(D, 2, "omitnan");
T.ShuffleRatio = T.diff2./T.ShuffleMedianDiff2;
T.EmpP = (1 + atLeast)./(1 + nScored);
T.ShuffleSigFrac = sum(sig, 2)./max(nScored, 1);
T = sortrows(T, ["EmpP", "diff2"], ["ascend", "descend"]);

medianRho = median(rho, "omitnan");
medianJac = median(jac, "omitnan");
insensitive = medianRho >= 0.9 && medianJac >= 0.8;
nEmp = nnz(T.EmpP <= 0.05);
% How many would get there by chance. EmpP takes the values k/(N+1), and
% under the null each is equally likely, so P(EmpP <= 0.05) is
% floor(0.05*(N+1))/(N+1): 1/21 at N = 20. Over hundreds of candidates that
% is a dozen or more pairs, so the count means something only beside it.
nChance = nCand*floor(0.05*(N + 1))/(N + 1);
% Not quoted in the verdict: random splits of one sample, where nothing
% differs, put two to four times this many pairs there (see the help).
caveat = [' EmpP is not calibrated for a change between conditions: random ' ...
    'splits of one sample give more pairs at EmpP <= 0.05 than an ' ...
    'exchangeable null predicts, so compare this count with ' ...
    'ten.xct2splitcalibration, and use q_value for calls.'];
% Two empty call sets "agree" completely; say so rather than quote a Jaccard
if ~any(sigReal) && ~any(sig, "all")
    calls = 'no pair is called at q < 0.05 in the real run or any shuffle';
else
    calls = sprintf('the q < 0.05 calls overlap %.0f%% (median Jaccard)', 100*medianJac);
end
if insensitive
    verdict = sprintf(['NETWORK-INSENSITIVE - with the networks'' gene labels ' ...
        'shuffled (the same way in both samples), the between-sample changes ' ...
        'still correlate %.3f (median Spearman of diff2) with the real run, and ' ...
        '%s. This differential result does not depend on the networks; do not ' ...
        'report it as network-derived. %d of %d pairs beat their shuffles at ' ...
        'EmpP <= 0.05.%s'], medianRho, calls, nEmp, nCand, caveat);
else
    verdict = sprintf(['NETWORK-SENSITIVE - shuffling the networks moves the ' ...
        'between-sample changes (median Spearman of diff2 %.3f with the real ' ...
        'run); %s. %d of %d pairs beat their shuffles at EmpP <= 0.05.%s'], ...
        medianRho, calls, nEmp, nCand, caveat);
end
if N < 19
    verdict = [verdict, sprintf([' (With %d shuffles EmpP cannot go below ' ...
        '%.3f; use NumShuffles >= 19 to read it at 0.05.)'], N, 1/(N + 1))];
end
if opts.Verbose
    fprintf('[xct2null] %s\n', verdict);
end

S = struct("Verdict", string(verdict), "NetworkSensitive", ~insensitive, ...
    "NumShuffles", N, "NumCandidates", nCand, "MedianRho", medianRho, ...
    "MedianSigJaccard", medianJac, "NumSigReal", nnz(sigReal), ...
    "NumEmpPBelow05", nEmp, "NumExpectedByChance", nChance, ...
    "MinEmpP", 1/(N + 1), "Rho", rho, "SigJaccard", jac);
end


function args = i_grnargs(grns, pS, pT)
% The four networks, CELLTYPE1's permuted by pS in both samples and
% CELLTYPE2's by pT, marked as already processed
args = {'grn_s1', grns.A_s1(pS, pS), 'grn_s1_processed', true, ...
    'grn_s2', grns.A_s2(pS, pS), 'grn_s2_processed', true, ...
    'grn_t1', grns.A_t1(pT, pT), 'grn_t1_processed', true, ...
    'grn_t2', grns.A_t2(pT, pT), 'grn_t2_processed', true};
end


function [T, grns] = i_xct2(sce1, sce2, a, b, args) %#ok<INUSD> used inside evalc
% One direction, database candidates, quietly
[~, T, grns] = evalc("ten.sctenifoldxct2(sce1, sce2, a, b, false, args{:})");
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
