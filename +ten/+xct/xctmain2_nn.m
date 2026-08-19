function [T] = xctmain2_nn(X_s1, X_t1, X_s2, X_t2, g, varargin)
% XCTMAIN2_NN  Differential cell-cell interaction, Path A (neural network).
%   Wraps ten.xct.xctmain_nn run on each sample independently, then performs
%   a differential fold-change test.  Requires Deep Learning Toolbox (R2022b+).
%
%   T = ten.xct.xctmain2_nn(X_s1, X_t1, X_s2, X_t2, g)
%   T = ten.xct.xctmain2_nn(X_s1, X_t1, X_s2, X_t2, g, Name, Value)
%
%   Inputs:
%     X_s1, X_t1 - genes × cells expression for sample 1 (source, target)
%     X_s2, X_t2 - genes × cells expression for sample 2 (source, target)
%     g          - gene name list matching rows (string or cell)
%
%   Name-Value pairs:
%     'twosided'  - run both directions (default: true)
%     'n_dim'     - embedding dimensions (default: 3, as the reference)
%     'mu'        - cross-type correspondence weight (default: 0.9, as merge.py)
%     'ncomp'     - PCNet principal components (default: 5)
%     'grn_q'     - GRN edge-filtering quantile; 0 disables (default: 0)
%     'n_steps'   - Adam training iterations (default: 1000)
%     'lr'        - Adam learning rate (default: 0.01)
%     'verbose'   - print progress (default: true)
%
%   Output:
%     T  - table: ligand, receptor, dist_s1, dist_s2, diff2, FC, p_value,
%          q_value, sorted by diff2 descending.
%          When twosided=true, T is a cell {T_dir1, T_dir2}.
%
%   THE TEST CHANGED ON 2026-07-20 AND SO DID `FC`. It now follows
%   stat.py::chi2_diff_test via TEN.I_XCTCHI2: diff2 = (dist_s1-dist_s2)^2 over
%   all gene pairs, scaled by its mean, right-tail chi-square, Benjamini-
%   Hochberg across all pairs. `FC` is the scaled statistic, NOT the former
%   dist_s1/dist_s2 ratio, so it is not comparable to values recorded before
%   that date, and direction is no longer encoded - read dist_s1 against
%   dist_s2 for that.

% ── Toolbox guard ─────────────────────────────────────────────────────────
if ~(exist('dlarray', 'builtin') || exist('dlarray', 'file'))
    error('xctmain2_nn:noDLT', ...
        ['Deep Learning Toolbox (R2022b+) is required for Path A. ' ...
         'Use ten.xct.xctmain2 or ten.sctenifoldxct2 (Path B, spectral) ' ...
         'as alternatives that do not need it.']);
end

p = inputParser;
addOptional(p, 'twosided', true,  @islogical);
% Differential path: merge.py sets n_dim=3 and mu=0.9, where the single-pair
% class uses mu=1.0. Follow merge.py here.
addOptional(p, 'n_dim',    3,     @(x) isnumeric(x) && x > 0);
addOptional(p, 'mu',       0.9,   @(x) isnumeric(x) && x > 0);
addOptional(p, 'ncomp',    5,     @(x) isnumeric(x) && x >= 2);
addOptional(p, 'grn_q',    0,     @(x) isnumeric(x) && x >= 0 && x < 1);
addOptional(p, 'n_steps',  1000,  @(x) isnumeric(x) && x > 0);
addOptional(p, 'lr',       0.01,  @(x) isnumeric(x) && x > 0);
addOptional(p, 'verbose',  true,  @islogical);
addParameter(p, 'w12mode', "outer", @(x) ismember(string(x), ["outer","lr"]));
addParameter(p, 'alpha',   0.5,   @(x) isnumeric(x) && x >= 0 && x <= 1);
addParameter(p, 'grnoffset', 1.0, @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(p, 'useparallel', false, @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
parse(p, varargin{:});
opts = p.Results;

pass = {'n_dim', opts.n_dim, 'mu', opts.mu, 'ncomp', opts.ncomp, ...
        'grn_q', opts.grn_q, 'n_steps', opts.n_steps, 'lr', opts.lr, ...
        'pval', 1.0, 'verbose', opts.verbose, ...
        'w12mode', opts.w12mode, 'alpha', opts.alpha, ...
        'grnoffset', opts.grnoffset, 'useparallel', opts.useparallel};

% ── Direction 1 ───────────────────────────────────────────────────────────
if opts.verbose, fprintf('[xctmain2_nn] Sample 1, direction 1...\n'); end
[~, Ea1] = ten.xct.xctmain_nn(X_s1, X_t1, g, 'twosided', false, pass{:});
if opts.verbose, fprintf('[xctmain2_nn] Sample 2, direction 1...\n'); end
[~, Ea2] = ten.xct.xctmain_nn(X_s2, X_t2, g, 'twosided', false, pass{:});
T1  = i_diff(Ea1, Ea2, g, opts.verbose, 'xctmain2_nn');

if opts.twosided
    if opts.verbose, fprintf('[xctmain2_nn] Sample 1, direction 2...\n'); end
    [~, Eb1] = ten.xct.xctmain_nn(X_t1, X_s1, g, 'twosided', false, pass{:});
    if opts.verbose, fprintf('[xctmain2_nn] Sample 2, direction 2...\n'); end
    [~, Eb2] = ten.xct.xctmain_nn(X_t2, X_s2, g, 'twosided', false, pass{:});
    T2  = i_diff(Eb1, Eb2, g, opts.verbose, 'xctmain2_nn');
    T   = {T1, T2};
else
    T = T1;
end

end % main


%% ── LOCAL FUNCTION ───────────────────────────────────────────────────────

function Td = i_diff(E1, E2, g, verbose, tag)
% I_DIFF  Differential chi-square test on aligned distances.
%   Follows merge.py's nn_aligned_diff together with stat.py::chi2_diff_test:
%   the statistic is diff2 = (dist_s1 - dist_s2)^2 evaluated over EVERY gene
%   pair, scaled by its own mean, referred to a right-tail chi-square, and
%   Benjamini-Hochberg corrected across all pairs. Restriction to the candidate
%   L-R pairs happens last.
%
%   The embeddings are required rather than the two candidate-only tables,
%   because mean(diff2) is taken over all pairs and it sets every p-value.
%
%   Replaced a two-tailed empirical test on log2(dist_s1/dist_s2) against
%   randomly cross-paired distances, which matched nothing in the reference.
%   See TEN.I_XCTCHI2.

if isempty(E1.P_s) || isempty(E2.P_s)
    Td = table(); return;
end
if ~isequal(E1.li_idx, E2.li_idx) || ~isequal(E1.ri_idx, E2.ri_idx)
    error([tag ':pairMismatch'], ...
        ['The two samples matched different L-R pairs, so their embeddings ' ...
         'cannot be compared elementwise. Both must be run on one gene list.']);
end

ng = size(E1.P_s, 1);
D1 = pdist2(double(E1.P_s), double(E1.P_t));
D2 = pdist2(double(E2.P_s), double(E2.P_t));

diff2_all = (D1 - D2).^2;                       % merge.py:110
[p_all, q_all, FC_all] = ten.i_xctchi2(diff2_all(:), dof=1, tail="right");

candLin = sub2ind([ng, ng], E1.li_idx, E1.ri_idx);
g = string(g(:));

Td = table(g(E1.li_idx), g(E1.ri_idx), D1(candLin), D2(candLin), ...
    diff2_all(candLin), FC_all(candLin), p_all(candLin), q_all(candLin), ...
    'VariableNames', {'ligand', 'receptor', 'dist_s1', 'dist_s2', 'diff2', ...
    'FC', 'p_value', 'q_value'});
Td = sortrows(Td, 'diff2', 'descend');          % stat.py:89

if verbose
    fprintf('[%s]   %d L-R pairs; %d with q < 0.05.\n', ...
        tag, height(Td), sum(Td.q_value < 0.05));
end
end
