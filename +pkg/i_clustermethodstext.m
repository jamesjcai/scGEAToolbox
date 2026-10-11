function [paragraph, refkeys] = i_clustermethodstext(method, options)
%I_CLUSTERMETHODSTEXT Methods text for cell clustering.
%
%   [paragraph, refkeys] = pkg.i_clustermethodstext(method, Name=Value)
%
%   Describes what SingleCellExperiment/clustercells did with METHOD, the
%   method ids it accepts ("louvainpc", "kmeans", "snndpc", "sc3", ...).
%
%   NAME-VALUE ARGUMENTS
%     NumClusters  clusters found
%     K            the requested number of clusters (embedding and
%                  expression methods; for louvainpc, the target the
%                  resolution was tuned to, or NaN)
%     Resolution   the Louvain resolution used (louvainpc)
%     NumPCs       principal components clustered (louvainpc)
%     Harmony      the PCs were the Harmony-corrected ones (louvainpc)
%     NumExtraGenes  genes added to the HVGs for the PCs (louvainpc)
%     NumDims      dimensions of the embedding clustered
%
%   See also SINGLECELLEXPERIMENT/CLUSTERCELLS, SC_LOUVAIN, PKG.I_METHODSREPORT.

arguments
    method (1,1) string
    options.NumClusters (1,1) double = NaN
    options.K (1,1) double = NaN
    options.Resolution (1,1) double = NaN
    options.NumPCs (1,1) double = 30
    options.Harmony (1,1) logical = false
    options.NumExtraGenes (1,1) double = 0
    options.NumDims (1,1) double = NaN
end

refkeys = "cai2020";
if isnan(options.NumDims)
    embedding = "the cell embedding";
else
    embedding = sprintf("the %d-dimensional cell embedding", options.NumDims);
end
k = sprintf("k = %d", options.K);

switch method
    case "louvainpc"
        if options.Harmony
            space = sprintf("the %d Harmony-corrected principal components", ...
                options.NumPCs);
        else
            space = sprintf("the top %d principal components of the 2,000 " + ...
                "most variable genes", options.NumPCs);
            if options.NumExtraGenes > 0
                space = space + sprintf(" and %d further selected genes", ...
                    options.NumExtraGenes);
            end
            space = space + " (library-size-normalized, log(x + 1)-" + ...
                "transformed and scaled per gene, clipped at 10)";
        end
        if isnan(options.K)
            resolution = sprintf("at resolution %g", options.Resolution);
        else
            resolution = sprintf("at resolution %.3g, tuned to give %d " + ...
                "clusters", options.Resolution, options.K);
        end
        paragraph = "Cells were clustered with the Louvain algorithm " + ...
            "(Blondel et al., 2008) on a shared-nearest-neighbour graph " + ...
            "(20 nearest neighbours by Euclidean distance; edges weighted by " + ...
            "the Jaccard overlap of the two neighbourhoods and pruned below " + ...
            "1/15) built from " + space + ", maximizing modularity with a " + ...
            "resolution parameter (Reichardt and Bornholdt, 2006) " + resolution;
        refkeys = [refkeys, "blondel2008", "reichardt2006"];
    case "louvain"
        paragraph = "Cells were clustered with the Louvain algorithm " + ...
            "(Blondel et al., 2008) on a shared-nearest-neighbour graph of " + ...
            embedding + ", with the resolution parameter (Reichardt and " + ...
            "Bornholdt, 2006) tuned to give " + sprintf("%d clusters", options.K);
        refkeys = [refkeys, "blondel2008", "reichardt2006"];
    case "kmeans"
        paragraph = "Cells were clustered by k-means (" + k + ", best of 5 " + ...
            "replicates, up to 1,000 iterations) on " + embedding;
    case "kmedoids"
        paragraph = "Cells were clustered by k-medoids (" + k + ") on " + embedding;
    case "spectclust"
        paragraph = "Cells were clustered by spectral clustering (MATLAB's " + ...
            "spectralcluster, " + k + ") on " + embedding;
    case "mbkmeans"
        paragraph = "Cells were clustered by mini-batch k-means with " + ...
            "k-means++ seeding (" + k + ", 5 replicates) on " + embedding;
    case "snndpc"
        paragraph = "Cells were clustered by shared-nearest-neighbour " + ...
            "density-peak clustering (SNN-DPC; Liu et al., 2018; " + k + ...
            ") on " + embedding;
        refkeys = [refkeys, "liu2018"];
    case "sc3"
        paragraph = "Cells were clustered by SC3 consensus clustering " + ...
            "(Kiselev et al., 2017), as implemented in scGEAToolbox, on the " + ...
            "2,000 most variable genes (" + k + ")";
        refkeys = [refkeys, "kiselev2017"];
    otherwise
        paragraph = "Cells were clustered with " + upper(method) + ...
            " on the 2,000 most variable genes (" + k + ")";
end
hitTarget = method == "louvainpc" && options.K == options.NumClusters;
if ~isnan(options.NumClusters) && ~hitTarget
    paragraph = paragraph + ", giving " + pkg.i_thousands(options.NumClusters) + " clusters.";
else
    paragraph = paragraph + ".";
end
end
