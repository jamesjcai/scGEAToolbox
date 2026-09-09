function [Y] = sc_ifft(X, num_CCs_to_modify)
%scGFT — Synthetic Cell Generation
% https://github.com/Sanofi-Public/PMCB-scGFT/blob/master/R/Core.R

% https://pubmed.ncbi.nlm.nih.gov/39843603/
% https://pubmed.ncbi.nlm.nih.gov/31919373/
% https://pubmed.ncbi.nlm.nih.gov/37699885/
% https://pubmed.ncbi.nlm.nih.gov/35510186/

if nargin < 2, num_CCs_to_modify = 10; end

[G, numc] = size(X);
% G = number of genes (also frequency components)
% numc = number of cells

wasSparse = issparse(X);
if wasSparse, X = full(X); end

% XRAW is kept because the back-conversion at the end needs the counts it
% started from. X used to be overwritten in place by the normalisation
% below and then handed to CONVERT_TO_RAW_COUNTS as its RAW_ORIG argument,
% so the "raw counts" it scaled were CP10K values. Measured on Poisson
% counts with library sizes of 107 to 428, the returned matrix had library
% sizes of 7940 to 10963 -- the CP10K scale, 20 to 90 times the input --
% and +gui/callback_GenerateSyntheticCells merges that straight into a real
% raw-count matrix, so the synthetic cells looked far deeper than the real
% ones. The commented-out cast at the bottom of this function confirms
% count-typed output was intended.
Xraw = X;

% A cell with no counts at all cannot be synthesised from: norm_libsize
% turns its library size into NaN, so its whole column becomes NaN.
emptyCell = sum(X, 1) == 0;
if any(emptyCell)
    warning('sc_ifft:emptyCells', ...
        ['%d of %d cell(s) have no counts. Their synthetic counterparts ', ...
        'are returned as all zeros and they are left out of the ', ...
        'component-variance estimate.'], nnz(emptyCell), numc);
end

X = pkg.norm_libsize(X, 1e4);
Xn = log1p(X);

% DFT along genes dimension (dim=1), per cell; ft_mtx is G-by-numc
ft_mtx = fft(Xn, [], 1);
assert(size(ft_mtx, 2) == numc)

a = abs(ft_mtx);            % amplitudes, G-by-numc

% Equation 20 is Xhat_j[k] = |X_j[k]| / || |X_j[:]| ||_2 -- an L2 norm over
% the FREQUENCIES of one cell. FT_MTX is G-by-numc with frequencies down
% the rows, so that norm runs along dimension 1. This was vecnorm(a, 2, 2),
% which normalises each frequency across cells instead: the columns then
% had norms of 1.0177 to 1.8029 rather than 1, and the per-component
% variance fed to A_k came out inflated by a median factor of 16.8, so
% every synthetic cell was perturbed far more violently than the method
% specifies. SCGFT_SYNTHESIZE in this repo implements the same equation on
% a cells-by-genes matrix, where the frequency axis is dimension 2, and
% takes its norm along that axis (lines 56-60); with the dimension
% corrected the two agree to 4.4e-16.
cellNorm = vecnorm(a, 2, 1);
cellNorm(cellNorm == 0) = 1;    % a cell with no signal keeps its zeros
a = a./cellNorm;

% OMITNAN because an empty cell's column is NaN, and VAR would otherwise
% return NaN for every frequency, which propagates through A_k into every
% cell's spectrum. That is how one empty cell in the selection used to turn
% all 30 synthetic cells in a 30-cell fixture into all zeros -- and
% invisibly, since max(0, NaN) is 0 in MATLAB rather than NaN.
sigma = var(a, 0, 2, 'omitnan');   % variance per frequency across cells

synthesized_data = zeros(G, numc);

valid_k = 1:floor((G-1)/2);
num_to_mod = min(num_CCs_to_modify, length(valid_k));

% There is no PARENT_IDX any more. It used to be
% randi([1 numc], numc, 1), so the cell used to convert a synthetic cell
% back to counts was picked at random and had nothing to do with the cell
% whose spectrum had been perturbed: on 30 cells only 3 were converted
% against their own parent and 13 original cells were never used as a
% parent at all. Column CN below is synthesised from column CN of FT_MTX,
% so column CN of XRAW is its parent, and CONVERT_TO_RAW_COUNTS now pairs
% them by position. SCGFT_SYNTHESIZE records meta.parent(row) = j, the same
% cell it perturbed.

for cn = 1:numc
    X_modified = ft_mtx(:, cn);     % G-by-1 column for this cell
    selected_k = valid_k(randperm(length(valid_k), num_to_mod));
    for k = selected_k
        if numc > 1
            A_k = 2 + sqrt(sigma(k+1)) * randn();
        else
            A_k = 2 + randn();
        end

        % Apply modification to k-th component and its conjugate pair
        % X'[k] = A_k * X[k] (Equation 23)
        % Note: MATLAB uses 1-based indexing, so k+1 for component k
        X_modified(k + 1) = A_k * X_modified(k + 1);

        % Conjugate pair at N-k (which is G-k in 1-based indexing)
        % X'[N-k] = A_k * X[N-k]
        conj_idx = G - k + 1;
        X_modified(conj_idx) = A_k * X_modified(conj_idx);
    end
    synthesized_data(:, cn) = real(ifft(X_modified));
end
synthesized_data = max(0, synthesized_data);

Y = convert_to_raw_counts(synthesized_data, Xraw, Xn);

% Hand the counts back in the caller's storage class. This block was
% commented out, so a sparse single count matrix came back as dense double:
% +gui/callback_GenerateSyntheticCells assigns the result into sub_sce.X and
% merges it with the rest of the dataset, so on a real gene panel that is a
% large dense block against a sparse one. XRAW rather than X, because X is
% the normalised matrix by this point -- which is what made the original
% block wrong as well as inert.
if wasSparse
    Y = sparse(cast(Y, 'like', Xraw));
else
    Y = cast(Y, 'like', Xraw);
end
end


function raw_synth = convert_to_raw_counts(synth_norm, raw_orig, norm_orig)
% Convert synthesized normalized data back to raw counts
% Uses relative change rate between synthesized and original normalized data
% All matrices are gene-by-cell (G x num_cells), and column I of each is the
% same cell: synthetic cell I is synthesised from original cell I, so they
% are paired by position rather than through a separate index.
%
% The conversion is exact when nothing was modified -- an untouched spectrum
% gives synth_norm == norm_orig, hence a relative change of 1 and the
% parent's counts back unchanged -- which is the invariant to hold on to.

    num_synth = size(synth_norm, 2);
    num_genes = size(synth_norm, 1);
    raw_synth = zeros(num_genes, num_synth);

    for i = 1:num_synth
        % Original normalized and raw values for this cell's parent, which
        % is this same column.
        orig_norm = norm_orig(:, i);   % G-by-1
        orig_raw  = raw_orig(:, i);    % G-by-1

        % Calculate relative change; avoid division by zero
        orig_norm_safe = orig_norm;
        orig_norm_safe(orig_norm_safe == 0) = 1e-10;

        relative_change = synth_norm(:, i) ./ orig_norm_safe;

        % Apply to raw counts
        raw_synth(:, i) = round(orig_raw .* relative_change);

        % Handle cases where original was zero
        zero_mask = (orig_norm == 0);
        raw_synth(zero_mask, i) = 0;
    end
    % Ensure non-negative integer counts
    raw_synth = max(0, round(raw_synth));
end


%{
https://chatgpt.com/s/t_699b270e1f2081919a22b558d6bb6e47

scGFT Synthetic Cell Generation

Input dataset: Lung_dev_scRNA (12,384 cells)

Generate:
  (•) 50% additional cells
  ( ) Custom number: [      ]

Mode:
  (•) Per cluster
  ( ) Global

Rare cell amplification:
  [x] Oversample clusters with < 2% cells

Output:
  (•) Append to dataset
  ( ) Create new dataset
%}
