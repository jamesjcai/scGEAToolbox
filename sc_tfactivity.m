function [cs, tflist, gcommon, numtargetgenes] = sc_tfactivity(X, g, ...
    Ttfgn, speciestag, methodid)
% SC_TFACTIVITY  Estimate transcription factor (TF) activity per cell.
%
% The activity level of a TF in a given cell is the extent to which it is
% exerting its regulatory potential on its target genes.
%
% USAGE:
%   [cs, tflist] = sc_tfactivity(X, g);
%   [cs, tflist] = sc_tfactivity(X, g, [], 'hs', methodid);
%
% INPUTS:
%   X         - genes-by-cells expression matrix (raw or normalized counts)
%   g         - gene name list (cell array or string array, length = rows of X)
%   Ttfgn     - (optional) custom TF-target table with columns {tf, target, mor};
%               if empty, the built-in DoRothEA database is used
%   speciestag - 'hs' (human, default) or 'mm' (mouse)
%   methodid  - scoring method (1–6, default = 1); see METHOD DETAILS below
%
% OUTPUTS:
%   cs             - m-by-n activity score matrix (m TFs × n cells)
%   tflist         - m-by-1 string array of TF names (same order as rows of cs)
%   gcommon        - genes used in scoring (intersection with database targets)
%   numtargetgenes - number of target genes used per TF
%
% -------------------------------------------------------------------------
% BACKGROUND
% -------------------------------------------------------------------------
% TF activity inference estimates how strongly a transcription factor is
% exerting its regulatory influence on its target genes in each cell —
% a protein-level quantity not directly observable from mRNA data alone.
%
% All six methods rely on the DoRothEA database of curated TF→target
% relationships, where each interaction is annotated with a signed
% "mode of regulation" (mor) weight:
%   +1  activating target (TF induces transcription)
%   -1  repressing target (TF suppresses transcription)
%   Fractional values encode regulatory confidence (levels A–E).
%
% Beyond this shared database, the methods diverge in two important ways:
%
%   REGULATORY SIGN HANDLING
%   ┌─────────────────────────────────────────────────────────────────┐
%   │ Methods 1–5:  use POSITIVE relationships only (mor > 0).        │
%   │   TF activity ≈ enrichment / expression of activated targets.   │
%   │   Repressing targets are discarded, losing half the information. │
%   │                                                                  │
%   │ Method 6 (VIPER):  uses the FULL signed regulon.                │
%   │   Both activating (mor > 0) and repressing (mor < 0) targets    │
%   │   contribute: high expression of activators OR low expression   │
%   │   of repressors both increase the TF activity score.            │
%   └─────────────────────────────────────────────────────────────────┘
%
%   AGGREGATION STRATEGY
%   Once the target gene set is defined, methods differ in how they
%   combine target gene expression values across cells into a scalar
%   activity score per TF per cell (see METHOD DETAILS below).
%
% Shared preprocessing steps (methods 1–5):
%   (1) Filter DoRothEA to mor > 0 entries; build TF-target matrix t.
%   (2) Intersect t columns with genes measured in X.
%   (3) Normalize X (log-library-size), except method 1 (rank-based).
%   (4) Aggregate target expression → activity score cs (nTF × nCells).
%
% VIPER (method 6) reloads the full table, builds a signed matrix, and
% uses a probit-rank enrichment statistic instead of simple aggregation.
%
% -------------------------------------------------------------------------
% METHOD DETAILS  (methodid 1–6, ordered fastest → slowest)
% -------------------------------------------------------------------------
% Method 1 — WMEAN (Weighted Mean)  [DEFAULT]
%   Fastest. Computes activity = t * X, where t is the TF-target weight
%   matrix (rows = TFs, cols = genes) and X is the normalized expression
%   matrix. Each score is the weighted sum of target gene expression per
%   cell, equivalent to a weighted mean when mor weights are uniform.
%   Cost: O(nTF · nGenes · nCells) — single matrix multiply.
%   Note: sensitive to library-size differences; normalization recommended.
%
%   NOT the same as Seurat's AddModuleScore (Tirosh et al., Science 2016),
%   which additionally subtracts a background score computed from
%   expression-bin-matched control genes. The implementation here omits
%   that correction step. In decoupleR nomenclature this is "WMEAN".
%
% Method 2 — UCell  (Andreatta & Carmona, Comput Struct Biotechnol J 2021)
%   PMID: 34285779  https://doi.org/10.1016/j.csbj.2021.06.043
%   Rank-based AUC scoring. For each cell, target genes are ranked by
%   expression. TF activity = AUC of target gene ranks relative to a
%   ceiling rank of 1500. Distribution-free and robust to outliers.
%   Formula: activity = 1 - U / (n × maxRank), where U is the rank-sum
%   statistic and n is the number of target genes.
%   Cost: O(nGenes · log(nGenes) · nCells) for ranking + O(nTF · nCells).
%   Note: rank-based; library-size normalization is NOT required.
%
% Method 3 — VIPER / aREA  (Alvarez et al., Nature Genetics 2016)
%   PMID: 27322546  https://doi.org/10.1038/ng.3593
%   Virtual Inference of Protein-activity by Enriched Regulon, using the
%   analytic Rank-based Enrichment Analysis (aREA) statistic.
%   Cost: O(nGenes · log(nGenes) · nCells) for ranking + vectorized NES.
%
%   KEY DIFFERENCE FROM ALL OTHER METHODS: uses the FULL signed regulon
%   (both activating and repressing targets with their mor weights).
%
%   Algorithm (per cell):
%     (a) Rank all G genes by expression (tiedrank, ascending).
%     (b) Map ranks to z-scores via the probit transform:
%           z_i = √2 · erfinv(2 · r_i/(G+1) − 1)  [≡ Φ^{-1}(r_i/(G+1))]
%         This converts the uniform rank distribution to N(0,1).
%     (c) Compute NES for TF t:
%           NES(t,c) = √n_t · Σ_i(w̃_i · z_i)
%         where w̃_i = mor_i / Σ|mor_j|  (L1-normalised signed weights),
%         n_t is the number of target genes with nonzero mor.
%   Interpretation: highly expressed activating targets (mor > 0) AND
%   lowly expressed repressing targets (mor < 0) both increase NES.
%   The √n_t factor stabilises variance across regulons of different sizes.
%   Note: rank-based; library-size normalization is NOT required.
%
% Method 4 — ULM / decoupleR  (Badia-i-Mompel et al., Bioinformatics 2022)
%   doi: 10.1093/bioinformatics/btac832
%   Univariate Linear Model. For each TF and each cell, fits a weighted
%   linear model: gene_expression ~ activity. Returns both activity scores
%   and p-values. The ONLY method providing per-cell statistical
%   significance. TFs with fewer than 5 target genes are auto-excluded.
%   Cost: O(nTF · nCells · nTargets) — loop over TFs and cells.
%
% Method 5 — Pseudo-inverse  (least-squares regression)
%   Solves activity = pinv(t') * X, finding the activity vector that best
%   reconstructs target gene expression via least squares. Accounts for
%   shared target genes across TFs.
%   Cost: O(nGenes · nTF²) for the SVD — expensive when nTF is large
%   (e.g. ~500 TFs and ~20 000 genes ≈ 5 × 10⁹ floating-point ops).
%
% Method 6 — NMF  (Cortal et al., Bioinformatics 2021)  🐢
%   PMID: 33135076  https://doi.org/10.1093/bioinformatics/btaa947
%   Slowest. Non-negative Matrix Factorization with the TF-target matrix W
%   fixed. Factorizes X ≈ W × H where W encodes regulon structure and H
%   (rows = TFs, cols = cells) is the activity matrix, optimized over 100
%   iterations using the beta-divergence (IS divergence, β = 0).
%   Cost: O(100 · nGenes · nTF · nCells) — iterative multiplicative updates.
%   Note: enforces non-negativity; most biologically interpretable.
%
% -------------------------------------------------------------------------
% REFERENCE DATABASE
% -------------------------------------------------------------------------
%   DoRothEA  (Garcia-Alonso et al., Genome Research 2019, PMID: 31340985)
%   https://doi.org/10.1101/gr.240663.118
%   Curated TF–target interactions with confidence levels A–E derived from
%   ChIP-seq, motif analysis, co-expression and literature curation.
%   Human: assets/DoRothEA_TF_Target_DB/dorothea_hs.mat
%   Mouse:  assets/DoRothEA_TF_Target_DB/dorothea_mm.mat

if nargin < 2, error('USAGE: [cs, tflist] = sc_tfactivity(X, g);'); end
if nargin < 5 || isempty(methodid), methodid = 1; end
if nargin < 4 || isempty(speciestag), speciestag = 'hs'; end
% if nargin < 3 || isempty(Ttfgn) % tf-by-gene matrix T from database
%     %folder=fileparts(mfilename('fullpath'));
%     %wrkpth=fullfile(folder,'assets',filesep,'DoRothEA_TF_Target_DB',filesep);
%     pw1 = fileparts(mfilename('fullpath'));
%     switch lower(speciestag)
%         case {'hs', 'human'}
%             %fname=[wrkpth 'dorothea_hs.mat'];
%             fname = fullfile(pw1, 'assets', 'DoRothEA_TF_Target_DB', 'dorothea_hs.mat');
%         case {'mm', 'mouse'}
%             %fname=[wrkpth 'dorothea_mm.mat'];
%             fname = fullfile(pw1, 'assets', 'DoRothEA_TF_Target_DB', 'dorothea_mm.mat');
%         otherwise
%             error('TF database is not available for the species.');
%
%     end
%         fprintf('\nReading ... %s.\n', fname);
%             load(fname, 'T');
%             Ttfgn = T(T.mor > 0, :);
%             fprintf('Only positive regulatory relationships are used.\n');
% end

if nargin < 3 || isempty(Ttfgn)
    pw1 = fileparts(mfilename('fullpath'));
    if strcmpi(speciestag, 'hs') || strcmpi(speciestag, 'human')
        fname = fullfile(pw1, 'assets', 'DoRothEA_TF_Target_DB', 'dorothea_hs.mat');
    elseif strcmpi(speciestag, 'mm') || strcmpi(speciestag, 'mouse')
        fname = fullfile(pw1, 'assets', 'DoRothEA_TF_Target_DB', 'dorothea_mm.mat');
    else
        error('TF database is not available for the species.');
    end

    if exist(fname, 'file')
        fprintf('\nReading ... %s.\n', fname);
        load(fname, 'T');
        if ~exist('T', 'var')
            error('Variable T not found in the file %s.', fname);
        end
    else
        error('File %s does not exist.', fname);
    end
    Ttfgn = T(T.mor > 0, :); % Filter positive regulatory relationships
    fprintf('Only positive regulatory relationships are used.\n');
end

if ~isnumeric(X) || ~ismatrix(X)
    error('X must be a numeric matrix.');
end
if ~iscellstr(g) && ~isstring(g)
    error('g must be a cell array of strings.');
end
% if ~all(isfield(Ttfgn.Properties.VariableNames, {'tf', 'target', 'mor'}))
if ~all(contains({'tf', 'target', 'mor'}, Ttfgn.Properties.VariableNames))
    error('Ttfgn must contain the fields: tf, target, and mor.');
end

try
    if issparse(X), X = full(X); end
catch
    warning('Keep using sparse X.');
end

if ~ismember(methodid, [2, 3]) % UCell (2) and VIPER (3) are rank-based; no normalization needed
    [X] = sc_norm(X);
    [X] = log1p(X);
end

% Save original X before it is trimmed to the positive-only regulon genes.
% VIPER (method 3) reloads the full signed regulon and intersects independently,
% so it must use X_orig rather than the trimmed X.
X_orig = X;

[gid, gnlist] = findgroups(string(Ttfgn.target));
[tid, tflist] = findgroups(string(Ttfgn.tf));
t = zeros(max(tid), max(gid));
t(sub2ind([max(tid), max(gid)], tid, gid)) = Ttfgn.mor;

fprintf('Using the Dorothea dadtabase that contains %d TFs and %d targets.\n', ...
    size(t, 1), size(t, 2));

[~, k, l] = intersect(upper(g), upper(gnlist));
t = t(:, l); % tf-by-gene matrix
X = X(k, :); % gene-by-cell matrix (trimmed to positive-only regulon genes)
fprintf('Using %d target genes that are present in the data.\n', size(t, 2));

if nargout > 2, gcommon = g(k); end

switch methodid
    case 1 % WMEAN — weighted mean of target gene expression
        cs = t * X;
        numtargetgenes = sum(t > 0, 2);

    case 2 % UCell — rank-based AUC  (see also: sc_cellscore)
        cs = zeros(size(t, 1), size(X, 2));
        R = tiedrank(-X);
        R(R > 1500) = 1500 + 1;
        numtargetgenes = zeros(size(t, 1), 1);
        for k = 1:size(t, 1)
            idx1 = t(k, :) > 0;
            n1 = sum(idx1);
            if n1 > 0
                u = sum(R(idx1, :)) - (n1 * (n1 - 1)) / 2;
                cs(k, :) = 1 - u / (n1 * 1500);
                numtargetgenes(k) = n1;
            end
        end
        cs(cs < 0) = 0;

    case 3 % VIPER / aREA — signed regulon, probit-rank enrichment
        disp('ref: Alvarez et al., Nature Genetics 2016.');
        disp('PMID: 27322546  https://doi.org/10.1038/ng.3593');
        disp('Using FULL signed regulon (activating AND repressing targets).');

        % ---- reload full regulon with both positive and negative mor ----
        if nargin < 3 || isempty(Ttfgn)
            pw1 = fileparts(mfilename('fullpath'));
            if strcmpi(speciestag, 'hs') || strcmpi(speciestag, 'human')
                fname = fullfile(pw1, 'assets', 'DoRothEA_TF_Target_DB', ...
                    'dorothea_hs.mat');
            elseif strcmpi(speciestag, 'mm') || strcmpi(speciestag, 'mouse')
                fname = fullfile(pw1, 'assets', 'DoRothEA_TF_Target_DB', ...
                    'dorothea_mm.mat');
            else
                error('TF database not available for species ''%s''.', speciestag);
            end
            load(fname, 'T');           % load full T (all mor signs)
            T_signed = T;
        else
            T_signed = Ttfgn;           % use provided table as-is (both signs)
        end

        % ---- build signed TF-target matrix (nTF × nGene) ---------------
        [gid2, gnlist2] = findgroups(string(T_signed.target));
        [tid2, tflist2] = findgroups(string(T_signed.tf));
        t2 = zeros(max(tid2), max(gid2));
        t2(sub2ind([max(tid2), max(gid2)], tid2, gid2)) = T_signed.mor;

        % ---- intersect with expression data -----------------------------
        % Use X_orig: the full (untrimmed) expression matrix, because the
        % shared preprocessing above already trimmed X to positive-only
        % regulon genes, which is a smaller set than the full signed regulon.
        [~, k2, l2] = intersect(upper(g), upper(gnlist2));
        t2 = t2(:, l2);               % nTF × nCommonGenes (signed)
        X2 = X_orig(k2, :);           % nCommonGenes × nCells (from full X)
        fprintf('VIPER using %d common target genes.\n', numel(k2));

        if nargout > 2, gcommon = g(k2); end

        G2 = size(X2, 1);  % number of genes after intersection
        N2 = size(X2, 2);  % number of cells

        % ---- (a) rank genes per cell (ascending: rank 1 = lowest expr) --
        R2 = zeros(G2, N2);
        for ci = 1:N2
            R2(:, ci) = tiedrank(X2(:, ci));
        end

        % ---- (b) probit transform: ranks → N(0,1) z-scores via erfinv ---
        %  Φ^{-1}(p) = √2 · erfinv(2p-1);  p = r/(G+1) avoids 0 and 1
        P2 = R2 / (G2 + 1);
        Z2 = sqrt(2) * erfinv(2 * P2 - 1);   % G2 × N2

        % ---- (c) compute NES for each TF --------------------------------
        nTF2 = size(t2, 1);
        cs   = zeros(nTF2, N2);
        numtargetgenes = zeros(nTF2, 1);
        for ki = 1:nTF2
            idx2 = t2(ki, :) ~= 0;
            n_t  = sum(idx2);
            if n_t > 0
                w      = t2(ki, idx2)';        % signed mor weights
                w_norm = w / sum(abs(w));       % L1-normalise within regulon
                % NES = √n · Σ(w̃_i · z_i)
                cs(ki, :) = sqrt(n_t) * sum(w_norm .* Z2(idx2, :), 1);
                numtargetgenes(ki) = n_t;
            end
        end
        tflist = tflist2;

    case 4 % ULM — univariate linear model (decoupleR)
        disp('ref: decoupleR: ensemble of methods that infer biological activities');
        disp('from omics data within a unified framework, Bioinformatics 2022.');
        disp('https://doi.org/10.1093/bioinformatics/btac832');
        % Build regulon struct from TF-target matrix t (nTF x nGenes)
        regulons_struct = struct('name', {}, 'targets', {}, 'weights', {});
        for ki = 1:size(t, 1)
            idx = find(t(ki, :) ~= 0);
            if ~isempty(idx)
                entry.name    = char(tflist(ki));
                entry.targets = idx(:);
                entry.weights = t(ki, idx)';
                regulons_struct(end+1) = entry; %#ok<AGROW>
            end
        end
        [ulm_activities, ulm_stats] = pkg.e_ULM(X, regulons_struct, ...
            'verbose', false);
        cs     = ulm_activities.scores;
        tflist = string(ulm_activities.regulatorNames)';
        numtargetgenes = ulm_stats.nTargetsUsed;

    case 5 % Pseudo-inverse — least-squares regression
        cs = pinv(t') * X;
        numtargetgenes = sum(t > 0, 2);

    case 6 % NMF — beta-divergence non-negative matrix factorization
        disp('ref: Bioinformatics, Volume 37, Issue 9, 1 May 2021, Pages 1234-1245,');
        disp('https://doi.org/10.1093/bioinformatics/btaa947');
        disp('PMID: 33135076 PMCID: PMC8189679 DOI: 10.1093/bioinformatics/btaa947');
        n = size(t, 1);
        v.WRfixed = n;
        v.W = t.';
        [~, cs] = NMF(X, n, 100, 0, v);
        numtargetgenes = sum(t > 0, 2);
end
end


function [W, H, bDsave] = NMF(V, R, Niter, beta, initialV)
% [W,H, bDsave] = NMF(V,R,Niter,beta,initialV)
%    NMF with beta divergence cost function.
% Input :
%   - V : power spectrogram to factorize (a MxN matrix)
%   - R : number of templates
%   - Niter : number of iterations
%   - beta (optional): beta used for beta-divergence (default : beta = 0, IS divergence)
%   - initialV (optional) : initial values of W, H (a struct with
%   fields W and H)
% Output :
%   - W : frequency templates (MxR array)
%   - H : temporal activation
%   - bDsave : evolution of beta divergence
%
% Copyright (C) 2010 Romain Hennequin

% https://github.com/romi1502/NMF-matlab
disp('NMF implementation: https://github.com/romi1502/NMF-matlab');

verbose = false;

eta = 1;

% size of input spectrogram
M = size(V, 1);
N = size(V, 2);

% initialization
if nargin == 5
    if isfield(initialV, 'H')
        H = initialV.H;
    else
        H = rand(R, N);
    end
    if isfield(initialV, 'W')
        W = initialV.W;
    else
        W = rand(M, R);
    end

    if isfield(initialV, 'HRfixed')
        HRfixed = initialV.HRfixed;
    else
        HRfixed = 0;
    end

    if isfield(initialV, 'WRfixed')
        WRfixed = initialV.WRfixed;
    else
        WRfixed = 0;
    end

else
    H = rand(R, N);
    W = rand(M, R);
    HRfixed = 0;
    WRfixed = 0;

    if nargin == 3
        beta = 0;
    end
end

% array to save the value of the beta-divergence
bDsave = zeros(Niter, 1);

% computation of Lambda (estimate of V) and of filters repsonse
Lambda = W * H;

% Waitbar
message = ['Computing NMF .... iteration : 0/', int2str(Niter), ' completed'];
h = waitbar(0, message);

% iterative computation
for iter = 1:Niter

    %     % compute beta divergence and plot its evolution
    bDsave(iter) = betaDiv(V + eps, Lambda + eps, beta);

    % update of W
    if not(WRfixed)
        W = W .* ((Lambda.^(beta - 2) .* V) * H' + eps) ./ ((Lambda.^(beta - 1)) * H' + eps);
    else
        W(:, WRfixed+1:end) = W(:, WRfixed+1:end) .* ((Lambda.^(beta - 2) .* V) * H(WRfixed+1:end, :)' + eps) ./ ((Lambda.^(beta - 1)) * H(WRfixed+1:end, :)' + eps);
    end

    % recomputation of Lambda (estimate of V)
    Lambda = W * H + eps;

    % update of H
    if not(HRfixed)
        H = H .* (W' * (Lambda.^(beta - 2) .* V) + eps) ./ (W' * (Lambda.^(beta - 1)) + eps);
    else
        H(1:HRfixed, :) = H(1:HRfixed, :) .* (W(:, 1:HRfixed)' * (Lambda.^(beta - 2) .* V) + eps) ./ (W(:, 1:HRfixed)' * (Lambda.^(beta - 1)) + eps);
    end
    % recomputation of Lambda (estimate of V)
    Lambda = W * H + eps;

    message = ['Computing NMF. iteration : ', int2str(iter), '/', int2str(Niter)];
    if verbose
        disp(message);
    end
    waitbar(iter / Niter, h, message);
end

% % normalization
% for r0=1:R
%     % normalization of templates
%     chosenNorm = 2;
%     normW = norm(W(:,r0),chosenNorm);
%     H(r0,:) = normW*H(r0,:);
%     W(:,r0) = W(:,r0)/normW;
% end

close(h)
% close
end


function bD = betaDiv(V, Vh, beta)
if beta == 0
    bD = sum((V(:) ./ Vh(:)) - log(V(:) ./ Vh(:)) - 1);
elseif beta == 1
    bD = sum(V(:) .* (log(V(:)) - log(Vh(:))) + Vh(:) - V(:));
else
    bD = sum(max(1 / (beta * (beta - 1)) * (V(:).^beta + (beta - 1) * Vh(:).^beta - beta * V(:) .* Vh(:).^(beta - 1)), 0));
end
end
