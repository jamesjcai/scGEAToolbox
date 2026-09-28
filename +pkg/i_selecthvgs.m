function keep = i_selecthvgs(X, g, numhvg)
%I_SELECTHVGS  Flag the NUMHVG most variable genes of a raw count matrix.
%
%   keep = pkg.i_selecthvgs(X, g, numhvg) returns a logical vector, one
%   entry per row of X, marking the genes to keep. X is raw counts (genes x
%   cells) and G the matching gene names. This is the gene set
%   SingleCellExperiment.embedcells feeds the embedders, and the one
%   PKG.E_CELLPCS takes principal components of.
%
%   See also SINGLECELLEXPERIMENT/EMBEDCELLS, PKG.E_CELLPCS, SC_ANALYTICFIT.

% Rank every gene most to least variable, then keep the top
% NUMHVG. SC_ANALYTICFIT scores each gene against the closed-form
% gamma-Poisson curve implied by the library sizes; SC_SPLINEFIT
% fits a smoothing spline through the gene cloud instead, and
% SC_HVG is the Brennecke CV^2 model. The analytic curve is
% defined at every mean, so unlike the spline it has no
% end-of-fit region where genes are scored by extrapolation.
%
% Scored on the labelled example sets by the share of each cell's
% 15 UMAP neighbours carrying its own annotated cell type, it
% beats SC_SPLINEFIT on a tight gene budget and ties at the 2000
% default: 0.957 against 0.947 (pancreas) and 0.962 against 0.952
% (leukocyte) at 300 genes; 0.967 against 0.969 and 0.967 against
% 0.966 at 2000. The two agree on only ~73% of the top 2000, so
% which one runs is not a detail.
%
% The other two stay as fallbacks in the order they were the
% default, newest first.
try
    Tranked = sc_analyticfit(X, g);
    granked = Tranked.genes;
catch ME
    warning(ME.message);
    try
        [~, ~, granked] = sc_splinefit(X, g, true, false, true);
    catch ME2
        warning(ME2.message);
        % NORMIT false so this ranks on the same raw counts the
        % other two rank on. SC_HVG defaults it to true, which
        % used to make an SC_SPLINEFIT failure switch scheme
        % silently, and to the worse of the two: 0.953 against
        % 0.924 at 2000 genes on the pancreas set, 0.837 against
        % 0.727 at 300, a tie on the leukocyte set. Library size
        % tracks cell identity wherever RNA content differs by
        % cell type, and normalising here discards that. NORMIT
        % changes the ranking only - SC_HVG returns raw counts
        % either way. Its gamma GLM may hit the iteration limit
        % on raw counts and warn; those scores were measured with
        % that warning firing.
        [~, ~, granked] = sc_hvg(X, g, true, false, ...
            false, false, true);
    end
end

% Top NUMHVG of what came back, not of what went in: all three
% rankers drop genes that are zero in every cell, so on a subset
% of the cells - one cell type isolated for subtype annotation,
% say - fewer genes come back than G holds.
%
% ISMEMBER against G rather than the ranker's own row order,
% because SC_ANALYTICFIT returns names only. Gene order does not
% matter downstream - PKG.I_APPENDGENES appends at the bottom
% either way, and an embedder treats genes as unordered features.
% A gene name appearing twice in G keeps both rows, so the
% count can exceed NUMHVG by the number of such duplicates; two
% rows under one name cannot be told apart here anyway.
nkeep = min(numhvg, numel(granked));
keep = ismember(g, granked(1:nkeep));
end
