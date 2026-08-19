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
parse(p, varargin{:});

cfg = ten.i_xct2cfg(p.Results, "sctenifoldxct2");

T = ten.i_xct2core(sce1, sce2, celltype1, celltype2, twosided, cfg);

end % sctenifoldxct2
