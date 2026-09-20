function [P, poolname, halfid] = i_poolpseudobulk(X, gid, coords)
%I_POOLPSEUDOBULK Summed raw counts per cell pool, and per half of each pool.
%
%   [P, poolname] = pkg.i_poolpseudobulk(X, gid, coords)
%   [P, poolname, halfid] = pkg.i_poolpseudobulk(___)
%
%   Builds the profiles that let a reference classifier label a pool without
%   labelling its cells. P has 3*npools columns in three blocks:
%
%       1:K        one profile per pool
%       K+1:2K     the same pools, first half only
%       2K+1:3K    the same pools, second half only
%
%   with K = numel(poolname) and the pools in sorted order throughout, so
%   column k, K+k and 2K+k all describe poolname(k).
%
%   Counts are summed, not averaged, and left un-normalised. That is what
%   SCimilarity's utils.get_centroid does - sum the raw counts, then
%   normalise that one profile to 1e4 and log1p it - and the normalisation
%   half of it belongs downstream, after the gene space has been aligned to
%   the model. Summing weights each cell by its sequencing depth rather than
%   equally, which is the vendor's choice, not an oversight here.
%
%   THE HALVES ARE THE POINT, not a bonus. A pseudobulk of a pool holding
%   two cell types is a profile in its own right: it can land confidently on
%   some third type, and nothing about that one profile reveals the problem.
%   Splitting the pool along its own first principal direction is what
%   separates two populations when there are two, so the halves disagreeing
%   with the whole is the signal that the pool should not be trusted. When
%   the pool is homogeneous the split is arbitrary and both halves reproduce
%   the whole. This is the cheap local analogue of SCimilarity's
%   query_coherence.
%
%   INPUTS:
%     X      - genes-by-cells raw counts
%     gid    - NumCells-by-1 pool ids, e.g. from PKG.I_PARTITIONCELLS
%     coords - NumCells-by-D coordinates the halves are split in, the same
%              space the pools were cut in
%
%   OUTPUTS:
%     P        - genes-by-3K summed counts, in the block layout above
%     poolname - K-by-1 string of pool ids, sorted
%     halfid   - NumCells-by-1, 1 or 2. A pool of 3 cells or fewer cannot
%                be bisected and is all half 1, so its two half profiles
%                are identical to each other and to the whole pool, and it
%                passes a halves-agree test by construction.
%
% See also PKG.I_PARTITIONCELLS, PKG.I_POOLEDANNOTATE.

arguments
    X {mustBeNonempty}
    gid {mustBeNonempty}
    coords (:,:) {mustBeNumeric, mustBeReal}
end

groups = string(gid(:));
ncells = numel(groups);
if size(X, 2) ~= ncells
    error('pkg:i_poolpseudobulk:sizeMismatch', ...
        'X has %d columns and GID has %d entries. Pass one pool id per cell.', ...
        size(X, 2), ncells);
end
if size(coords, 1) ~= ncells
    error('pkg:i_poolpseudobulk:coordsMismatch', ...
        'COORDS has %d rows and GID has %d entries. Pass one row per cell.', ...
        size(coords, 1), ncells);
end

[poolname, ~, poolindex] = unique(groups);
npools = numel(poolname);

% Cells are visited by pool through one sort, not with FIND(poolindex==k)
% per pool: at one pool per 50 cells that test is quadratic in the cell
% count, and would cost more than the classifier this is meant to spare.
[~, order] = sort(poolindex);
poolsize = accumarray(poolindex, 1);
last = cumsum(poolsize);
first = last - poolsize + 1;

halfid = ones(ncells, 1);
for k = 1:npools
    rows = order(first(k):last(k));
    if numel(rows) < 4
        continue;   % cannot be bisected; see HALFID in the help above
    end
    % TargetSize = half the pool makes I_PARTITIONCELLS stop after exactly
    % one split: it splits a node larger than 1.5*TargetSize and each half
    % is below that, so two groups come back for every pool of 4 or more.
    sub = pkg.i_partitioncells(coords(rows, :), ...
        TargetSize=ceil(numel(rows)/2));
    if max(sub) ~= 2
        error('pkg:i_poolpseudobulk:notBisected', ...
            'Pool "%s" of %d cells split into %d parts, not 2.', ...
            poolname(k), numel(rows), max(sub));
    end
    halfid(rows) = sub;
end

% One sparse matrix multiply per block rather than a loop over pools.
P = [i_sum(X, poolindex, npools, true(ncells, 1)), ...
    i_sum(X, poolindex, npools, halfid == 1), ...
    i_sum(X, poolindex, npools, halfid == 2)];
end


function S = i_sum(X, poolindex, npools, keep)
% Sum the columns of X within each pool, over the cells KEEP selects.

ncells = numel(poolindex);
selector = sparse(1:ncells, poolindex, double(keep), ncells, npools);
if isa(X, 'single')
    % X is single sparse from R2025a on, and a double selector would make
    % the product double - or refuse to multiply at all.
    selector = single(selector);
end
S = X*selector;
end
