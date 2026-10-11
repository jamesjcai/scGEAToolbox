function [c_cell_type_tx, cL_cell_type_tx, T] = sc_celltypeanno(X, g, c_cluster_id, species, options)
%SC_CELLTYPEANNO Marker-based cell type annotation of clusters (PanglaoDB).
%
%   [c_cell_type_tx, cL_cell_type_tx, T] = SC_CELLTYPEANNO(X, g, c_cluster_id)
%   [...] = SC_CELLTYPEANNO(X, g, c_cluster_id, species)
%   [...] = SC_CELLTYPEANNO(..., Name=Value)
%
%   Labels each cluster with the PanglaoDB cell type whose marker genes it
%   expresses most strongly, by cell-type activity (CTA) scoring
%   (Franzen et al. 2019, Database baz046; RUN.ML_ALONA). Deterministic:
%   the same input always gives the same labels.
%
%   This is one annotation method among several. To choose a method -
%   marker-based (this one), LLM-based, reference-based, or a consensus of
%   several - and have the labels and their provenance written to a
%   SingleCellExperiment, use SC_ANNOTATECELLS, whose default
%   Method="markers" calls this function.
%
%   INPUTS
%     X             genes-by-cells count matrix. Raw counts are expected;
%                   it is library-size normalized and log2(x + 1)
%                   transformed internally.
%     g             gene names, one per row of X, as a string array. Matched
%                   to PanglaoDB case-insensitively.
%     c_cluster_id  cluster index per cell, integers 1..K
%     species       "human" (default) or "mouse": which PanglaoDB sheet of
%                   assets/PanglaoDB/celltypes.xlsx supplies the markers
%
%   NAME-VALUE ARGUMENTS
%     Summary  "median" (default, as the paper) or "mean"; see below
%     Alpha    FDR above which a cluster is "Unknown" (default 0.05)
%     Verbose  print a Methods paragraph describing the run, with
%              references, to paste into a manuscript (default true)
%
%   OUTPUTS
%     c_cell_type_tx  - one label per cell
%     cL_cell_type_tx - one label per cluster, for clusters 1..K
%     T               - one row; columns Ck_Cell_Type and Ck_CTA_Score for
%                       each cluster k
%
%   SCORING
%   For each cell type, the median log-normalized expression of each of its
%   markers in the cluster, weighted so that markers shared by many cell
%   types count less, is summed and divided by the cube root of the number
%   of markers found. The highest-scoring type wins, unless a one-sided
%   Fisher test finds its markers not over-represented among the genes the
%   cluster expresses (BH-adjusted p > 0.05): the cluster is then labelled
%   "Unknown". See RUN.ML_ALONA for the details.
%
%   The median follows the paper and has a cost on droplet data with
%   strong ambient RNA: a transcript present in every cell (islet hormones
%   in a pancreas run) has a high median in every cluster, while a real
%   marker detected in half the cells has a median of zero. On
%   example_data/new_example_sce.mat this turns one fibroblast cluster into
%   "Enteroendocrine cells". Summary="mean" scores with the mean instead;
%   that is a deviation from the paper, and the printed Methods text says
%   so.
%
%   EXAMPLE
%     [X, g] = sc_read10xdir('path/to/10x_data');
%     cid = sc_cluster_s(sc_umap(X, 2), 10);
%     [ct, cL, T] = sc_celltypeanno(X, g, cid, "mouse");
%
%   See also SC_ANNOTATECELLS, RUN.ML_ALONA, LLM.E_CELLTYPEANNO,
%   SC_CELLTYPEANNOREF.
arguments
    X
    g
    c_cluster_id
    species = 'human'
    options.Summary (1,1) string {mustBeMember(options.Summary, ["median", "mean"])} = "median"
    options.Alpha (1,1) double {mustBePositive} = 0.05
    options.Verbose (1,1) logical = true
end
if isempty(species), species = 'human'; end
[T] = run.ml_alona(X, g, c_cluster_id, 'species', species, 'bestonly', true, ...
    'summary', options.Summary, 'alpha', options.Alpha);
cL_cell_type_tx = table2array(T(:, 1:2:end))';
c_cell_type_tx = cL_cell_type_tx(c_cluster_id);
if options.Verbose
    pkg.i_ctamethodstext(Species=string(species), Summary=options.Summary, ...
        Alpha=options.Alpha, NumClusters=numel(cL_cell_type_tx), ...
        NumUnknown=sum(string(cL_cell_type_tx) == "Unknown"));
end
end
