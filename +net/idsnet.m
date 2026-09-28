function A = idsnet(X)
% Construct co-expression network using the InterDependence Score (IDS)
%
% A = net.idsnet(X)
%
% X - genes x cells, library-size normalised and log1p transformed
%     (what gui.i_transformx returns); not raw UMI counts
% A - genes x genes IDS matrix, values in [0, 1], diagonal set to 0
%
% IDS detects nonlinear dependence, including dependence confined to a small
% subpopulation of cells, at roughly the cost of a correlation matrix.
%
% Each gene is min-max scaled into [0, 1] before scoring. The kernel has a
% fixed bandwidth centred at zero, and on log-normalised values spanning
% 0 to ~5 the high-order Taylor features are dominated by a few extreme
% cells. Simulated UMI data (2000 cells, 100 genes, a 20-gene program on in
% 5% of cells, 3-4 fold depth spread) gave AUROC 0.50 without scaling,
% 0.92 with [0, 1], 0.79 with [0, 2], and 0.91 for Pearson; raw UMI counts
% gave 0.66, their null inflated by shared sequencing depth.
%
% Scores use PNorm=2 (root mean square of the 36 feature correlations)
% rather than the reference default "max": 0.94 vs 0.92 in the 5% case
% above and 0.71 vs 0.68 with the program on in 1% of cells.
% ref: Radhakrishnan et al., PNAS 2025. DOI:10.1073/pnas.2509860122
%
% See also: pkg.e_ids, net.minet, net.xicornet, net.distcorrnet

arguments
    X {mustBeNumeric}
end

A = pkg.e_ids(X.', ScaleRange=[0 1], PNorm=2);
A(1:size(A, 1) + 1:end) = 0;
end
