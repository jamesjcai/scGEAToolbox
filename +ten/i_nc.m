function [XM] = i_nc(X, nsubsmpl, ncom, csubsmpl, usebootstrp, useparallel, parseed)
% NC - network construction
%
% input: X -  n (genes/features) x m (cells/samples) matrix
% output XM - k multi-layer network array (n x n x k)
%
% Networks are built by TEN.I_XCTGRN, shared with the scTenifoldXct entry
% points, with symmetrize=false because the tensor decomposition downstream
% needs them directed.
%
% EDGE FILTERING WAS BROKEN AND IS NOW FIXED (2026-07-20). This loop read:
%
%     XM(:, :, k) = ten.e_filtadjc(A, 0.95, false);
%     a = max(abs(A(:)));
%     XM(:, :, k) = A ./ a;
%
% The second assignment overwrote the first, so the 0.95 filter was computed on
% every subsample and then thrown away. The two lines were competing
% treatments, not a sequence: e_filtadjc normalizes by the top-10 mean AND
% filters, while the surviving line normalized by the max and did not filter.
%
% The reference settles it. scTenifoldNet's makeNetworks defaults to q = 0.95,
% and its pcNet scales by max(abs(.)) and then zeroes entries below the qth
% quantile, without symmetrizing. So filtering was intended and the surviving
% line was the incomplete half. Networks are now built with q = 0.95.
%
% THIS CHANGES RESULTS. The networks go from fully dense to roughly 5% of their
% edges, which is the point, but any scTenifoldNet or scTenifoldKnk output
% recorded before this date was computed on unfiltered networks.
%
% One deviation remains, and it is small and fully accounted for.
% ten.e_filtadjc takes its quantile over the NONZERO entries, where the
% reference takes one over all entries. The only zeros in a dense pcNet output
% are the n diagonal ones, so e_filtadjc keeps 5% of n^2-n rather than of n^2:
% measured on a 600-gene network, 17970 edges against the reference's 18000,
% with e_filtadjc's set a strict subset of it (Jaccard 0.9983). The 30-edge
% gap is exactly 0.05*600.
%
% e_filtadjc also rescales by the top-10 mean rather than the max. That is a
% uniform factor, so it does not change WHICH edges survive - only their
% magnitudes, by a per-layer constant.
%
% REPRODUCIBILITY (fixed 2026-07-20). This function used to call
% rng("shuffle") before the subsampling loop, which did two unwanted things:
% it discarded whatever seed the caller had set, so scTenifoldNet and
% scTenifoldKnk could not be reproduced even from a fixed seed, and it left the
% global stream reseeded on exit, altering every later random draw in the
% session. Neither is appropriate in a library function; the toolbox's other
% computational users of the RNG, such as ten.i_nnembed, preserve the caller's
% stream instead.
%
% The call is gone. The caller's stream now governs, so
%
%     rng(42); XM = ten.i_nc(X);
%
% is repeatable, while consecutive calls without reseeding still draw different
% subsamples because the stream advances normally. The only remaining effect on
% the RNG is that consumption which any sampling function has.
%
% Note this makes runs repeatable from a fresh MATLAB session, which always
% starts from the same default seed - previously each session differed.
%
% PARALLEL SUBSAMPLING (added 2026-08-20). USEPARALLEL runs the subsample loop
% as a parfor, one worker per network. It is OFF by default and the serial
% branch below is byte-for-byte the code that was there before, so nothing
% changes for any existing caller.
%
% READ THIS BEFORE EXPECTING A 10x SPEEDUP. The serial loop is not
% single-threaded: TEN.I_XCTGRN's header records that multithreaded BLAS
% already gives the per-network SVD loop 3.0x on 600 genes and 5.7x on 1856
% (20 cores, 784 cells). Process-pool workers get ONE computational thread
% each, so a 10-worker parfor trades that BLAS parallelism for task
% parallelism rather than adding to it. On 20 cores the honest expectation is
% 10 cores of single-threaded work against roughly 6 cores' worth of BLAS
% speedup - a modest gain, not an order of magnitude. Measure it on the
% workload that matters before planning around it; the gain grows with
% nsubsmpl and shrinks as the gene count rises and BLAS scales better.
%
% RNG: THE PARALLEL PATH DRAWS DIFFERENT SUBSAMPLES. parfor cannot consume the
% global stream in loop order, so the parallel branch gives each network its
% own substream of a single PARSEED. Network k then depends only on
% (PARSEED, k) - reproducible from a seed, and independent of how many workers
% ran or in what order they finished, which is a stronger guarantee than the
% serial path has. But it is NOT the same draw the serial path makes: there,
% subsample k depends on everything the stream consumed before it, including
% the svds calls inside I_XCTGRN. The two branches are statistically
% equivalent and neither is more correct. They are not interchangeable
% mid-analysis: results from one should not be pooled with results from the
% other without noting it, exactly as two different seeds would not be.
%
% PARSEED defaults to a draw from the caller's stream, so rng(42) before the
% call still makes the whole thing reproducible, while consecutive unseeded
% calls still differ.
import ten.*

if nargin < 5, usebootstrp = true; end % using m-out-of-n bootstrap (false by default)
% using jackknife (by default)
if nargin < 4, csubsmpl = 500; end % number of cells in subsamples
if nargin < 3, ncom = 3; end % number of components for PC regression
if nargin < 2, nsubsmpl = 10; end % number of subsamples
if nargin < 6 || isempty(useparallel), useparallel = false; end
if nargin < 7 || isempty(parseed), parseed = randi(intmax('int32')); end

n  = size(X, 1);
n0 = size(X, 2);
XM = zeros(n, n, nsubsmpl);

% The bootstrap fallback is a property of the population, not of the
% iteration, so it is decided once. The serial loop below used to re-evaluate
% it every pass and latch it to true; same outcome, stated once.
if n0 < csubsmpl * 1.15
    usebootstrp = true;
end

if ~useparallel
    % ---- serial: unchanged ------------------------------------------------
    for k = 1:nsubsmpl
        fprintf('Building network...%d of %d\n', k, nsubsmpl);
        if usebootstrp % bootstrap
            i = randi(n0, 1, csubsmpl);
            Xrep = X(:, i);
        else % jackknife
            Xrep = X(:, randperm(n0));
            Xrep = Xrep(:, 1:csubsmpl);
        end
        % Shared PCNet builder, same code as the scTenifoldXct entry points
        % use. q=0.95 and symmetrize=false match makeNetworks' defaults in the
        % reference R package (nComp=3, scaleScores=TRUE, symmetric=FALSE,
        % q=0.95).
        XM(:, :, k) = ten.i_xctgrn(Xrep, ncom, 0.95, false, false, ...
            symmetrize=false);
    end
    return
end

% ---- parallel -------------------------------------------------------------
% Draw every subsample up front, each from its own substream, then slice the
% columns before the loop. Slicing here rather than indexing X inside the
% parfor matters: X as a broadcast variable would be copied whole to every
% worker, whereas csubsmpl columns per iteration is the only part any worker
% needs.
Xrep = cell(nsubsmpl, 1);
for k = 1:nsubsmpl
    st = RandStream.create('mrg32k3a', 'NumStreams', nsubsmpl, ...
        'StreamIndices', k, 'Seed', parseed);
    if usebootstrp
        Xrep{k} = X(:, randi(st, n0, 1, csubsmpl));
    else
        r = randperm(st, n0);
        Xrep{k} = X(:, r(1:csubsmpl));
    end
end

fprintf('Building %d networks in parallel (parseed=%d)...\n', nsubsmpl, parseed);
parfor k = 1:nsubsmpl
    XM(:, :, k) = ten.i_xctgrn(Xrep{k}, ncom, 0.95, false, false, ...
        symmetrize=false);
end
end
