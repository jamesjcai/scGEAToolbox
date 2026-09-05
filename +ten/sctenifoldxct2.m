function [T] = sctenifoldxct2(sce1, sce2, celltype1, celltype2, twosided, varargin)
% SCTENIFOLDXCT2  Differential cell-cell interaction across two samples.
%   Native MATLAB implementation (Path B, PCNet GRNs).  No Python required.
%
%   T = ten.sctenifoldxct2(SCE1, SCE2, CT1, CT2)
%   T = ten.sctenifoldxct2(SCE1, SCE2, CT1, CT2, TWOSIDED)
%   T = ten.sctenifoldxct2(SCE1, SCE2, CT1, CT2, TWOSIDED, Name, Value)
%
%   Identifies L-R pairs whose interaction strength differs significantly
%   between two conditions (SCE1 vs SCE2) for the same pair of cell types.
%   For each direction, spectral manifold alignment is run independently on
%   each sample; the squared difference of embedding distances is then referred
%   to a chi-square distribution, as merge.py's chi2_diff_test does.
%
%   Inputs:
%     SCE1, SCE2  - SingleCellExperiment objects for the two conditions
%     CT1         - source cell type label (string)
%     CT2         - target cell type label (string)
%     TWOSIDED    - logical, run CT1→CT2 AND CT2→CT1 (default: true)
%
%   Name-Value pairs:
%     'n_dim'   - spectral embedding dimensions (default: 3, as the reference)
%     'mu'      - cross-type correspondence weight (default: 0.9, as merge.py)
%     'verbose' - print progress messages (default: true)
%     'w12mode' - "lr" (default) builds the correspondence from cognate
%                 ligand-receptor pairs; "outer" reproduces the reference's
%                 dense expression outer product, in which the L-R database
%                 plays no part in the alignment. The same mode is applied to
%                 both samples, which is what keeps their distances
%                 comparable. See TEN.I_XCTW12.
%     'alpha'   - mean/variance blend for w12mode="outer" (default: 0.5)
%     'grnoffset' - constant added to every element of each GRN block before
%                 assembly (default: 1.0, as core.py). See TEN.I_XCTBLOCK.
%                 Until 2026-07-20 this function applied no offset at all,
%                 having kept its own inline block assembly; it now shares
%                 TEN.I_XCTBLOCK with the other entry points.
%     'useparallel' - build the GRNs with a parfor (default: false, and usually
%                 the slower choice; see TEN.I_XCTGRN)
%
%   Precomputed-GRN Name-Value pairs:
%     'grn_s1'   - ng-by-ng adjacency to use for CT1 in SCE1, instead of
%                  building one via net.pcrnet, where ng = numel(sce1.g).
%                  Must be in the same gene order as sce1.g (== sce2.g; the
%                  two are required to share one gene list). [] (default)
%                  builds as usual.
%     'grn_t1'   - same, for CT2 in SCE1.
%     'grn_s2'   - same, for CT1 in SCE2.
%     'grn_t2'   - same, for CT2 in SCE2.
%                  Any of the four may be supplied independently. Passing one
%                  skips that sample/cell-type's pcrnet call, but the matrix
%                  still goes through TEN.I_XCTGRN's scale/filter/symmetrize,
%                  so mu and grnoffset stay meaningful. Substituting a network
%                  built by a different method changes what this differential
%                  test is comparing, not just how it got there - see
%                  TEN.I_XCTGRN and TEN.SCTENIFOLDXCT's grn1/grn2, which this
%                  mirrors for the two-sample case.
%     'grn_s1_processed', 'grn_t1_processed', 'grn_s2_processed',
%     'grn_t2_processed' - true if the matching grn_* is itself a grns.A_s/
%                  A_t TEN.SCTENIFOLDXCT previously returned - skips scale/
%                  filter/symmetrize a second time, which is NOT a no-op for
%                  the q=0.75 filter (default false; see TEN.I_XCTGRN's
%                  'processed' option).
%
%   Discovery-mode Name-Value pairs:
%     'candidates' - "database" (default) restricts the output to L-R
%                  database matches, as published. "all" instead ranks
%                  EVERY gene pair by |diff2| and reports the top 'topN',
%                  regardless of database membership - for finding pairs
%                  whose difference between samples the database does not
%                  already know to look for. The chi-square test already
%                  runs over every gene pair internally either way (merge.py
%                  builds its full table before testing), so this changes
%                  only which rows are returned, plus adds an is_known_lr
%                  column. Combine with w12mode="outer" for a run that does
%                  not use the L-R database at all. See TEN.I_XCT2CORE's
%                  i_xct2discover.
%     'topN'     - how many top pairs to report when candidates="all"
%                  (default 200). Ignored when candidates="database".
%
%   Output:
%     T  - table: ligand, receptor, dist_s1, dist_s2, diff2, FC, p_value,
%          q_value, sorted by diff2 descending.
%            diff2   : (dist_s1 - dist_s2)^2, the test statistic
%            FC      : dof*diff2/mean(diff2) over ALL gene pairs, the value
%                      referred to the chi-square
%            p_value : right-tail chi-square probability; small means the pair's
%                      distance changed more between conditions than typical
%            q_value : Benjamini-Hochberg across all gene pairs
%          When twosided=true, T is a cell array {T_ct1_ct2, T_ct2_ct1}.
%
%   THE TEST CHANGED ON 2026-07-20 AND SO DID `FC`. This function previously
%   computed FC = dist_s1/dist_s2 and a two-tailed empirical p-value for
%   |log2(FC)| against randomly cross-paired gene distances. That matched
%   nothing in the reference. It now follows stat.py::chi2_diff_test via
%   TEN.I_XCTCHI2. `FC` is therefore no longer a ratio and is not comparable to
%   values recorded before that date; a ratio near 1 used to mean "no change",
%   whereas a small FC now means the same thing. Direction is no longer encoded
%   at all, because the statistic is squared - read dist_s1 against dist_s2 for
%   that.
%
%   The engine moved to TEN.I_XCT2CORE on 2026-07-26 so that this function and
%   TEN.SCTENIFOLDXCT2_GLYCO cannot drift apart, mirroring what TEN.I_XCTCORE
%   does for the single-sample pair. Numerically nothing changed, with one
%   exception: SCE1 and SCE2 are now required to carry the same gene list in
%   the same order, which was always assumed and is now checked.
%
%   This is a MATLAB replacement for +run/py_scTenifoldXct2.
%   Reference: Ma et al., Cell Systems 2023. PMID:36787742
%
% see also: TEN.I_XCT2CORE, TEN.SCTENIFOLDXCT2_GLYCO, TEN.I_XCTCHI2,
%           TEN.SCTENIFOLDXCT

if nargin < 5, twosided = true; end

p = inputParser;
% n_dim follows merge.py; mu stays 0.9 because merge_scTenifoldXct sets 0.9
% where the single-pair class uses 1.0. The reference is not self-consistent
% here, and this is the differential path, so it follows merge.py.
addOptional(p, 'n_dim',   3,    @(x) isnumeric(x) && x > 0);
addOptional(p, 'mu',      0.9,  @(x) isnumeric(x) && x > 0);
addOptional(p, 'verbose', true, @islogical);
addParameter(p, 'w12mode', "lr",  @(x) ismember(string(x), ["lr","outer"]));
addParameter(p, 'alpha',   0.5,   @(x) isnumeric(x) && x >= 0 && x <= 1);
addParameter(p, 'grnoffset', 1.0, @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(p, 'useparallel', false, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
addParameter(p, 'grn_s1', [], @(x) isempty(x) || (isnumeric(x) && ismatrix(x)));
addParameter(p, 'grn_t1', [], @(x) isempty(x) || (isnumeric(x) && ismatrix(x)));
addParameter(p, 'grn_s2', [], @(x) isempty(x) || (isnumeric(x) && ismatrix(x)));
addParameter(p, 'grn_t2', [], @(x) isempty(x) || (isnumeric(x) && ismatrix(x)));
addParameter(p, 'grn_s1_processed', false, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
addParameter(p, 'grn_t1_processed', false, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
addParameter(p, 'grn_s2_processed', false, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
addParameter(p, 'grn_t2_processed', false, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
addParameter(p, 'candidates', "database", @(x) ismember(string(x), ["database","all"]));
addParameter(p, 'topN', 200, @(x) isnumeric(x) && isscalar(x) && x > 0);
parse(p, varargin{:});

cfg = ten.i_xct2cfg(p.Results, "sctenifoldxct2");
cfg.candidates = string(p.Results.candidates);
cfg.topN = round(p.Results.topN);

ng = numel(sce1.g);
i_checkgrnsize(p.Results.grn_s1, ng, 'grn_s1');
i_checkgrnsize(p.Results.grn_t1, ng, 'grn_t1');
i_checkgrnsize(p.Results.grn_s2, ng, 'grn_s2');
i_checkgrnsize(p.Results.grn_t2, ng, 'grn_t2');
cfg.grn_s1 = double(p.Results.grn_s1);
cfg.grn_t1 = double(p.Results.grn_t1);
cfg.grn_s2 = double(p.Results.grn_s2);
cfg.grn_t2 = double(p.Results.grn_t2);
cfg.grn_s1_processed = logical(p.Results.grn_s1_processed);
cfg.grn_t1_processed = logical(p.Results.grn_t1_processed);
cfg.grn_s2_processed = logical(p.Results.grn_s2_processed);
cfg.grn_t2_processed = logical(p.Results.grn_t2_processed);

T = ten.i_xct2core(sce1, sce2, celltype1, celltype2, twosided, cfg);

end % sctenifoldxct2


%% ---- validate a precomputed GRN against the gene count ----
function i_checkgrnsize(A, ng, name)
if isempty(A), return; end
if ~isequal(size(A), [ng, ng])
    error("TEN:SCTENIFOLDXCT2:BadGrnSize", ...
        "'%s' must be %d-by-%d (sce1.g has %d genes). Received %s.", ...
        name, ng, ng, ng, mat2str(size(A)));
end
end % i_checkgrnsize
