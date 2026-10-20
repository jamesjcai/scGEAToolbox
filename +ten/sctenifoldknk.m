function [T, A0, genelist] = sctenifoldknk(X, genelist, kogene, varargin)
% T=sctenifoldknk(X,genelist,"Foxp3");
%
% X is a gene x cell matrix from the wild-type
%
% The knockout removes KOGENE's outgoing edges, as in R; see TEN.I_KNK for the
% 2026-09-25 fix of its direction. A0 is returned in NET.PCRNET's orientation
% (row k = regression of gene k), over the genes of the third output.
%
% THE TOOLBOX VARIANT IS THE DEFAULT ('reference', false): libsize + log1p,
% no QC, CP (ten.i_cpalstt) of rank NSUBSMPL, normalized Laplacian, and at
% least 10 outgoing links on the KO gene. The third output is then GENELIST.
%
% 'ncom' (default 20) is the number of principal components each pcNet
% regression in TEN.I_NC uses. It was hardcoded at 3 until 2026-10-19 and
% changed after pipelines/tenifold_laplacian (README "pcNet components"):
% on a linear-SEM simulation scored against interventional truth, the
% median responder AUROC was 0.64 at ncom = 3 and 0.91 at 20 with 25
% regulator modules, 0.90 and 1.00 with 5; ncom = 20 beat 3 in all twelve
% seed x regime cells. On three real populations it made knockouts more
% specific (top-50 overlap between different knockouts 0.27-0.46 -> 0.21-
% 0.27). Three components cannot represent more than a few co-regulated
% modules, so each gene's regression borrows from the same few axes. Network
% construction costs about 4x as much (1.7 s -> 6.5 s at 1000 genes, 10
% networks). Pass 'ncom', 3 to reproduce earlier results. This is a
% different bottleneck from the one in TEN.KNKNULLCONTROL's RANK COLLAPSE
% note: the CP rank still caps A0 at rank NSUBSMPL whatever ncom is.
% It is capped at one below the gene count. Reference mode ignores 'ncom'
% (R uses 3).
%
% LITE MODE ('lite', true) builds A0 as llm.run_sctenifoldknk does: libsize
% + log1p, then one TEN.I_PCNET on all cells with q = 0, unsymmetrized and
% unscaled - no subsampling, no edge filter, no CP. 'nsubsmpl', 'csubsmpl',
% 'smplmethod' and 'tdmethod' are ignored; 'ncom' applies. It is opt-in, but
% on the evidence of pipelines/tenifold_laplacian (README "Which step loses
% the signal") it is the better network: in the linear-SEM simulation it
% matched |Pearson| (median responder AUROC 0.87 / 0.72 at weak signal, base
% 0.52 / 0.56), on three real populations it gave the most knockout-specific
% output (top-50 overlap between knockouts 0.10-0.12, base 0.21-0.27), and it
% builds 5-10x faster. Subsampling to 500 cells is what loses the signal in
% simulation; the q = 0.95 filter plus CP is what makes different knockouts
% report the same hub genes on real data. A0 is dense, so every gene passes
% the 10-link check, and its rank is not capped by CP (see TEN.KNKNULLCONTROL).
%
%     [T, A0, g] = ten.sctenifoldknk(X, genelist, "Foxp3", 'lite', true);
%
% WHY THE TOOLBOX VARIANT IS THE DEFAULT. Reference mode was the default for
% part of 2026-09-25 and was switched back after a paired benchmark
% (pipelines/tenifold_laplacian, README "Why reference mode scores lower"):
% on planted co-expression changes R's pipeline scored AUROC 0.72 against
% 0.93 for the toolbox variant, the loss coming from R's rounding of the denoised network to one
% decimal, its rank-3 decomposition and CPM without log1p. Use reference mode
% when the numbers have to agree with R or scTenifoldpy.
%
% REFERENCE MODE ('reference', true) reproduces scTenifoldKnk() of the R
% package (1.1.4) run with its defaults as scTenifoldKnk(X, gKO, seed = SEED):
%
%     [T, A0, g] = ten.sctenifoldknk(X, genelist, "Foxp3", 'reference', true);
%
% X must then be raw counts. Every step follows R: TEN.I_SCQC for QC and CPM
% (no log1p; KOGENE must pass QC), TEN.I_NCREF with q = 0.9, TEN.I_CPALS
% with K = 3 rounded to 3 decimals, and TEN.I_KNK in its reference mode
% (unnormalized, unsymmetrized alignment in 2 dimensions; TEN.I_DR reference
% output, sorted by p). 'seed' (default 1) is R's seed; 'smplmethod' and
% 'tdmethod' are ignored, and so is 'sorttable': the table is always sorted,
% as R's is. A0 has a zero diagonal. G, the genes after QC, is what A0 and T
% are over - fewer than GENELIST when QC drops any, so index A0 with G, not
% with GENELIST. Results match R to floating-point precision, except for
% edges tied at the quantile threshold or values right on a rounding
% boundary.
%
% CHECK THE RESULT AGAINST A NULL BEFORE USING IT. T's p-values come from
% TEN.I_DR, which fits a chi-square to drdist/mean(drdist). That asks whether
% a gene moved more than the average gene moved in this run - not whether it
% moved because of KOGENE. On a network with strong hub structure the same
% genes top the ranking whichever row is zeroed, and the enrichment computed
% from them describes the network rather than the knockout.
%
% TEN.KNKNULLCONTROL settles it by knocking out random genes from the same A0
% and rescoring against that background. It takes ~13 s at 3000 genes, needs the A0
% returned here, and reports a verdict:
%
%     [T, A0, g] = ten.sctenifoldknk(X, genelist, "Foxp3");
%     [Tnull, S] = ten.knknullcontrol(A0, "Foxp3", g);
%     disp(S.Verdict)
%
% see also: TEN.KNKNULLCONTROL, TEN.I_KNK, TEN.I_DR, TEN.SCTENIFOLDNET
import ten.*

if nargin < 3
    error(sprintf('USAGE: T=sctenifoldknk(X,genelist,kogene);\n       T=sctenifoldnet_m(X0,X1,genelist,''qqplot'',true);'));
end
if isscalar(kogene) && isnumeric(kogene)
    idx = kogene;
else
    idx = find(genelist == kogene, 1);
    if isempty(idx)
        error("KOGENE should be a member of GENELIST.");
    end
end
p = inputParser;
addOptional(p, 'sorttable', false, @islogical);
addOptional(p, 'smplmethod', "bootstrap", @(x) (isstring(x) | ischar(x)) & ismember(lower(string(x)), ["jackknife", "bootstrap"]));
addOptional(p, 'tdmethod', "CP", @(x) (isstring(x) | ischar(x)) & ismember(upper(string(x)), ["CP", "TUCKER"]));
addOptional(p, 'nsubsmpl', 10, @(x) fix(x) == x & x > 0);
addOptional(p, 'csubsmpl', 500, @(x) fix(x) == x & x > 0);
addOptional(p, 'ncom', 20, @(x) isscalar(x) && fix(x) == x && x > 0);
addOptional(p, 'savegrn', false, @islogical);
addOptional(p, 'reference', false, @islogical);
addOptional(p, 'lite', false, @islogical);
addOptional(p, 'seed', 1, @(x) isnumeric(x) && isscalar(x) && fix(x) == x);
parse(p, varargin{:});
if p.Results.reference && p.Results.lite
    error("'reference' and 'lite' build different networks. Set at most one of them to true.");
end
if p.Results.reference
    [T, A0, genelist] = i_reference(X, genelist, idx, p.Results);
    return
end
dosort = p.Results.sorttable;
tdmethod = p.Results.tdmethod;
nsubsmpl = p.Results.nsubsmpl;
csubsmpl = p.Results.csubsmpl;
% Capped below the gene count, which pcNet needs, so the default still runs
% on small gene sets
ncom = min(p.Results.ncom, size(X, 1) - 1);
smplmethod = p.Results.smplmethod;
savegrn = p.Results.savegrn;

switch upper(tdmethod)
    case "CP"
        tdmethod = 1;
    case "TUCKER"
        tdmethod = 2;
end
switch lower(smplmethod)
    case "jackknife"
        usebootstrp = false;
    case "bootstrap"
        usebootstrp = true;
end

if size(X, 1) ~= length(genelist)
    error('Length of genelist should be the same as the number of rows of X0 or X1.');
end

X = sc_norm(X, "type", "libsize");
X = log1p(X);

if p.Results.lite
    A0 = ten.i_pcnet(X, ncom, 0, false, false, symmetrize=false, scale=false);
else
    [XM] = i_nc(X, nsubsmpl, ncom, csubsmpl, usebootstrp);
    [A0] = i_td1(XM, tdmethod);
end
    if savegrn
        tstr = matlab.lang.makeValidName(strrep(sprintf("GRN_created_on_% s", datetime)," ", "_at_"));
        save(tstr, 'A0', 'genelist', '-v7.3', '-nocompression');
        fprintf('\nConstructed gene regulatory network (GRN) is saved in %s.mat\n', tstr);
    end
T = ten.i_knk(A0, idx, genelist, dosort);
end


function [T, A0, genelist] = i_reference(X, genelist, idx, opts)
% scTenifoldKnk() of the R package; see REFERENCE MODE in the header
if size(X, 1) ~= numel(genelist)
    error('Length of genelist should be the same as the number of rows of X.');
end
genelist = string(genelist(:));
kogene = genelist(idx);
[X, genelist] = ten.i_scqc(X, genelist);
idx = find(genelist == kogene, 1);
if isempty(idx)
    error('ten:sctenifoldknk:koGeneRemoved', ...
        'The KO gene %s is not present in the count matrix after quality control.', kogene);
end

XM = ten.i_ncref(X, 'NumNets', opts.nsubsmpl, 'NumCells', opts.csubsmpl, ...
    'NumComp', 3, 'Q', 0.9, 'Seed', opts.seed);
A0 = ten.i_cpals(XM, 3, NumDecimals=3, Seed=opts.seed);
clear XM
A0(1:(size(A0, 1) + 1):end) = 0;
if opts.savegrn
    tstr = matlab.lang.makeValidName(string(datetime));
    save(sprintf('A0_%s', tstr), 'A0', 'genelist', '-v7.3', '-nocompression');
end
T = ten.i_knk(A0, idx, genelist, true, 0, true);
end
