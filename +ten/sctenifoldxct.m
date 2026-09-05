function [T, grns] = sctenifoldxct(sce_ori, celltype1, celltype2, twosided, varargin)
% SCTENIFOLDXCT  Native MATLAB implementation of scTenifoldXct.
%   Detects ligand-receptor-mediated cell-cell interactions by spectral
%   manifold alignment of gene regulatory networks.  No Python required.
%
%   T = ten.sctenifoldxct(SCE, CELLTYPE1, CELLTYPE2)
%   T = ten.sctenifoldxct(SCE, CELLTYPE1, CELLTYPE2, TWOSIDED)
%   T = ten.sctenifoldxct(SCE, CELLTYPE1, CELLTYPE2, TWOSIDED, Name, Value)
%   [T, grns] = ten.sctenifoldxct(...)   % also return the GRNs actually used
%
%   Inputs:
%     SCE       - SingleCellExperiment object
%     CELLTYPE1 - source cell type label (string)
%     CELLTYPE2 - target cell type label (string)
%     TWOSIDED  - logical, run both directions (default: true)
%
%   Name-Value pairs:
%     'n_dim'   - spectral embedding dimensions (default: 3, as the reference)
%     'mu'      - cross-cell-type correspondence weight (default: 1.0, as the reference)
%     'pval'    - p-value threshold for null test, 1.0 returns all pairs
%                 (default: 1.0)
%     'verbose' - print progress messages (default: true)
%
%   Correspondence Name-Value pairs:
%     'w12mode' - "lr" (default) builds the correspondence block from cognate
%                 ligand-receptor pairs, weight 1 each. "outer" reproduces the
%                 reference implementation's dense expression outer product
%                 over all gene pairs, in which the L-R database plays no part
%                 in the alignment and is used only to name candidates. See
%                 TEN.I_XCTW12.
%     'alpha'   - mean/variance blend for w12mode="outer" (default: 0.5)
%     'grnoffset' - constant added to EVERY element of each GRN block before
%                 assembly (default: 1.0, as core.py). This is what connects
%                 the graph and makes the Laplacian PSD; it also makes the
%                 alignment graph dense, so subset genes for large problems.
%                 Pass 0 to reproduce results recorded before it was restored.
%                 See TEN.I_XCTBLOCK.
%
%   Solver Name-Value pairs:
%     'solver'   - "spectral" (default) or "nn". The spectral route solves the
%                  alignment in closed form via the Laplacian eigenvectors.
%                  The "nn" route trains two networks to minimise the same
%                  loss, as the published Python implementation does; because
%                  its embedding is confined to the image of that network
%                  family, it is a different method rather than a different
%                  solver, and generally attains a higher loss. Requires Deep
%                  Learning Toolbox.
%     'n_steps'  - Adam iterations for solver="nn" (default: 1000)
%     'lr'       - Adam learning rate for solver="nn" (default: 0.01)
%     'seed'     - RNG seed for solver="nn" (default: 0)
%
%   Precomputed-GRN Name-Value pairs:
%     'grn1'     - ng-by-ng adjacency to use for CELLTYPE1 instead of building
%                  one via net.pcrnet, where ng = numel(sce.g). Must be in the
%                  same gene order as sce.g. [] (default) builds as usual.
%     'grn2'     - same, for CELLTYPE2.
%                  Passing either skips that side's pcrnet call - the
%                  expensive step at cell counts in the tens of thousands -
%                  but the matrix still goes through TEN.I_XCTGRN's scale/
%                  filter/symmetrize so cfg.mu and cfg.grnoffset stay
%                  meaningful. A network built by a different method (a
%                  different ncomp, a bootstrapped/denoised one, ...) is a
%                  different quantity than the one net.pcrnet(X,3) would have
%                  built here, not a faster way to get the same answer - see
%                  TEN.I_XCTGRN.
%     'grn1_processed', 'grn2_processed' - true if the matching grn1/grn2 is
%                  itself a grns.A_s/A_t this function previously returned
%                  (see the grns output below) - skips scale/filter/symmetrize
%                  a second time, which is NOT a no-op for the q=0.75 filter
%                  (default false; see TEN.I_XCTGRN's 'processed' option).
%
%   Discovery-mode Name-Value pairs:
%     'candidates' - "database" (default) restricts the output to L-R
%                  database matches, as published. "all" instead ranks
%                  EVERY gene pair by embedding distance and reports the
%                  closest 'topN', regardless of whether either gene is a
%                  known ligand or receptor - for finding candidate pairs
%                  the database does not already know about. Combine with
%                  w12mode="outer" for a run that does not use the L-R
%                  database at all; w12mode="lr" (the default) still uses it
%                  to build the alignment even when candidates="all" only
%                  changes what gets reported afterward. The output schema
%                  differs in this mode: no p_value (there is no separate
%                  background left once every pair is a candidate) - instead
%                  a percentile column (smaller = closer than more of the
%                  other ~25 million gene pairs) and is_known_lr, marking
%                  which of the top pairs happen to already be in the
%                  database. See TEN.I_XCTCORE's i_xctdiscover.
%     'topN'     - how many top pairs to report when candidates="all"
%                  (default 200). Ignored when candidates="database".
%
%   Output:
%     T    - candidates="database" (the default): table with columns ligand,
%            receptor, dist, correspondence, p_value. candidates="all":
%            table with columns ligand, receptor, dist, percentile,
%            is_known_lr. If twosided=true, T is a cell array {T1, T2} where
%            T1 is CELLTYPE1->CELLTYPE2 and T2 is CELLTYPE2->CELLTYPE1.
%     grns - requested by asking for a second output; struct with fields
%            A_s, A_t (the ng-by-ng adjacency actually used for CELLTYPE1/
%            CELLTYPE2 - freshly built, or grn1/grn2 passed through TEN.
%            I_XCTGRN's scale/filter/symmetrize either way) and genes. One
%            pair regardless of twosided - both directions align the same
%            two networks, just swapped - so there is nothing to return
%            twice. Building it costs nothing extra: A_s/A_t already exist
%            by the time T does. This is how to keep a network net.pcrnet
%            just built (cell counts in the tens of thousands make that
%            expensive) for later reuse - e.g. as another call's grn1/grn2 -
%            without rebuilding it.
%
%   Algorithm:
%     Builds partial-correlation GRNs for each cell type via net.pcrnet, then
%     performs spectral manifold alignment via graph Laplacian eigenvectors on
%     the combined [GRN_source, W_LR; W_LR', GRN_target] weight matrix.
%     Ligand-receptor pairs from the built-in database act as correspondences
%     that pull ligand and receptor gene embeddings together.  Significance is
%     assessed by a nonparametric left-tail null test comparing candidate L-R
%     distances to a sampled background.
%
%   For glycan-context-aware correspondences, see TEN.SCTENIFOLDXCT_GLYCO,
%   which shares this function's engine and adds the glyco channels without
%   altering anything here.
%
%   This function is a drop-in MATLAB replacement for +run/py_scTenifoldXct.
%   Reference: Ma et al., Cell Systems 2023. PMID:36787742
%
% see also: TEN.SCTENIFOLDXCT_GLYCO, TEN.I_XCTCORE, TEN.I_ALIGNEMBED

if nargin < 4, twosided = true; end

% The glyco options used to live here. They now have their own entry point, so
% point the caller at it rather than letting inputParser report an unrecognised
% name-value pair.
i_rejectglyco(varargin);

p = inputParser;
addOptional(p, 'n_dim',   3,    @(x) isnumeric(x) && x > 0);
addOptional(p, 'mu',      1.0,  @(x) isnumeric(x) && x > 0);
addOptional(p, 'pval',    1.0,  @(x) isnumeric(x) && x >= 0);
addOptional(p, 'verbose', true, @islogical);
addParameter(p, 'w12mode', "lr",       @(x) ismember(string(x), ["lr","outer"]));
addParameter(p, 'alpha',   0.5,        @(x) isnumeric(x) && x >= 0 && x <= 1);
addParameter(p, 'grnoffset', 1.0,      @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(p, 'useparallel', false, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
addParameter(p, 'solver',  "spectral", @(x) ismember(string(x), ["spectral","nn"]));
addParameter(p, 'n_steps', 1000,       @(x) isnumeric(x) && x > 0);
addParameter(p, 'lr',      0.01,       @(x) isnumeric(x) && x > 0);
addParameter(p, 'seed',    0,          @isnumeric);
addParameter(p, 'grn1',    [],         @(x) isempty(x) || (isnumeric(x) && ismatrix(x)));
addParameter(p, 'grn2',    [],         @(x) isempty(x) || (isnumeric(x) && ismatrix(x)));
addParameter(p, 'grn1_processed', false, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
addParameter(p, 'grn2_processed', false, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
addParameter(p, 'candidates', "database", @(x) ismember(string(x), ["database","all"]));
addParameter(p, 'topN', 200, @(x) isnumeric(x) && isscalar(x) && x > 0);
parse(p, varargin{:});

cfg = ten.i_xctcfg(p.Results, "sctenifoldxct");
cfg.candidates = string(p.Results.candidates);
cfg.topN = round(p.Results.topN);

% -- Subset to the two cell types -----------------------------------------
sce = copy(sce_ori);
idx = sce.c_cell_type_tx == celltype1 | sce.c_cell_type_tx == celltype2;
sce = sce.selectcells(idx);

ng = numel(sce.g);
i_checkgrnsize(p.Results.grn1, ng, 'grn1');
i_checkgrnsize(p.Results.grn2, ng, 'grn2');
cfg.grn1 = double(p.Results.grn1);
cfg.grn2 = double(p.Results.grn2);
cfg.grn1_processed = logical(p.Results.grn1_processed);
cfg.grn2_processed = logical(p.Results.grn2_processed);

X = sce.X;
if issparse(X), X = full(X); end
X = single(X);

if nargout > 1
    [T, grns] = ten.i_xctcore(X, sce.g, sce.c_cell_type_tx, celltype1, ...
        celltype2, twosided, cfg);
else
    T = ten.i_xctcore(X, sce.g, sce.c_cell_type_tx, celltype1, celltype2, ...
        twosided, cfg);
end

end % sctenifoldxct


%% ---- validate a precomputed GRN against the gene count ----
function i_checkgrnsize(A, ng, name)
if isempty(A), return; end
if ~isequal(size(A), [ng, ng])
    error("TEN:SCTENIFOLDXCT:BadGrnSize", ...
        "'%s' must be %d-by-%d (sce.g has %d genes after subsetting to " + ...
        "the two cell types). Received %s.", name, ng, ng, ng, mat2str(size(A)));
end
end % i_checkgrnsize


%% ---- explanatory error for the relocated glyco options ----
function i_rejectglyco(args)
moved = ["glyco", "glycochannel", "glycolambda", "glycomode", "glycoscale"];
for k = 1:2:numel(args)
    if ~(ischar(args{k}) || isstring(args{k})), continue; end
    name = lower(string(args{k}));
    if ismember(name, moved)
        error("TEN:SCTENIFOLDXCT:GlycoMoved", ...
            "'%s' is no longer accepted here. The glyco channels moved to " + ...
            "TEN.SCTENIFOLDXCT_GLYCO, which takes the same arguments minus " + ...
            "'glyco' itself. Replace the call with " + ...
            "ten.sctenifoldxct_glyco(sce, celltype1, celltype2, ...).", name);
    end
end
end % i_rejectglyco
