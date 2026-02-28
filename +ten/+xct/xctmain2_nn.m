function [T] = xctmain2_nn(X_s1, X_t1, X_s2, X_t2, g, varargin)
%XCTMAIN2_NN  Differential cell-cell interaction, Path A (neural network).
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
%     'n_dim'     - embedding dimensions (default: 50)
%     'mu'        - cross-type correspondence weight (default: 0.9)
%     'corr_thr'  - |Pearson r| cut-off for GRN edges (default: 0.3)
%     'n_steps'   - Adam training iterations (default: 1000)
%     'lr'        - Adam learning rate (default: 0.01)
%     'verbose'   - print progress (default: true)
%
%   Output:
%     T  - table: ligand, receptor, dist_s1, dist_s2, FC, p_value
%          FC < 1: interaction stronger in sample 1
%          FC > 1: interaction stronger in sample 2
%          When twosided=true, T is a cell {T_dir1, T_dir2}.
%
%   The differential p-value uses the same empirical null as xctmain2:
%   cross-pairing d1 and d2 distances at random simulates the null of
%   exchangeable sample labels.

% ── Toolbox guard ─────────────────────────────────────────────────────────
if ~(exist('dlarray', 'builtin') || exist('dlarray', 'file'))
    error('xctmain2_nn:noDLT', ...
        ['Deep Learning Toolbox (R2022b+) is required for Path A. ' ...
         'Use ten.xct.xctmain2 (Path B Lite) or ten.sctenifoldxct2 ' ...
         '(Path B) as toolbox-free alternatives.']);
end

p = inputParser;
addOptional(p, 'twosided', true,  @islogical);
addOptional(p, 'n_dim',    50,    @(x) isnumeric(x) && x > 0);
addOptional(p, 'mu',       0.9,   @(x) isnumeric(x) && x > 0);
addOptional(p, 'corr_thr', 0.3,   @(x) isnumeric(x) && x >= 0 && x <= 1);
addOptional(p, 'n_steps',  1000,  @(x) isnumeric(x) && x > 0);
addOptional(p, 'lr',       0.01,  @(x) isnumeric(x) && x > 0);
addOptional(p, 'verbose',  true,  @islogical);
parse(p, varargin{:});
opts = p.Results;

pass = {'n_dim', opts.n_dim, 'mu', opts.mu, 'corr_thr', opts.corr_thr, ...
        'n_steps', opts.n_steps, 'lr', opts.lr, ...
        'pval', 1.0, 'verbose', opts.verbose};

% ── Direction 1 ───────────────────────────────────────────────────────────
if opts.verbose, fprintf('[xctmain2_nn] Sample 1, direction 1...\n'); end
Ta1 = ten.xct.xctmain_nn(X_s1, X_t1, g, 'twosided', false, pass{:});
if opts.verbose, fprintf('[xctmain2_nn] Sample 2, direction 1...\n'); end
Ta2 = ten.xct.xctmain_nn(X_s2, X_t2, g, 'twosided', false, pass{:});
T1  = i_diff(Ta1, Ta2, opts.verbose);

if opts.twosided
    if opts.verbose, fprintf('[xctmain2_nn] Sample 1, direction 2...\n'); end
    Tb1 = ten.xct.xctmain_nn(X_t1, X_s1, g, 'twosided', false, pass{:});
    if opts.verbose, fprintf('[xctmain2_nn] Sample 2, direction 2...\n'); end
    Tb2 = ten.xct.xctmain_nn(X_t2, X_s2, g, 'twosided', false, pass{:});
    T2  = i_diff(Tb1, Tb2, opts.verbose);
    T   = {T1, T2};
else
    T = T1;
end

end % main


%% ── LOCAL FUNCTION ───────────────────────────────────────────────────────

function Td = i_diff(T_s1, T_s2, verbose)
%I_DIFF  Join two single-sample XCT tables and compute differential FC test.

    if isempty(T_s1) || isempty(T_s2) || height(T_s1) == 0 || height(T_s2) == 0
        Td = table(); return;
    end

    % ── Join on (ligand, receptor) ────────────────────────────────────────
    key1 = upper(string(T_s1.ligand)) + "|" + upper(string(T_s1.receptor));
    key2 = upper(string(T_s2.ligand)) + "|" + upper(string(T_s2.receptor));
    [~, ia, ib] = intersect(key1, key2, 'stable');

    if isempty(ia)
        Td = table();
        if verbose
            warning('xctmain2_nn:noCommonPairs', ...
                'No common L-R pairs between the two samples.');
        end
        return;
    end

    lig = T_s1.ligand(ia);
    rec = T_s1.receptor(ia);
    d1  = T_s1.dist(ia);
    d2  = T_s2.dist(ib);

    % ── Fold change ───────────────────────────────────────────────────────
    FC  = d1 ./ max(d2, 1e-8);
    lfc = log2(FC);

    % ── Null: cross-pair random FC (H0: d1 and d2 exchangeable) ──────────
    n_pairs = numel(d1);
    n_null  = max(50000, 100 * n_pairs);
    rng_state = rng;
    ri1 = randi(n_pairs, n_null, 1);
    ri2 = randi(n_pairs, n_null, 1);
    rng(rng_state);
    null_lfc  = log2(d1(ri1) ./ max(d2(ri2), 1e-8));

    % Two-tailed p-value
    p_vals = sum(abs(null_lfc) >= abs(lfc)', 1)' ./ n_null;

    % ── Output table ─────────────────────────────────────────────────────
    Td = table(lig, rec, d1, d2, FC, p_vals, ...
        'VariableNames', {'ligand','receptor','dist_s1','dist_s2','FC','p_value'});
    Td = sortrows(Td, 'p_value', 'ascend');

    if verbose
        fprintf('[xctmain2_nn]   %d common L-R pairs; %d with p < 0.05.\n', ...
            height(Td), sum(Td.p_value < 0.05));
    end
end
