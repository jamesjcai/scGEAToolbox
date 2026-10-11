function [paragraph, refkeys] = i_embedmethodstext(part, options)
%I_EMBEDMETHODSTEXT Methods text for feature selection and embedding.
%
%   [paragraph, refkeys] = pkg.i_embedmethodstext("features", Name=Value)
%   [paragraph, refkeys] = pkg.i_embedmethodstext(method, Name=Value)
%
%   "features" describes the genes SingleCellExperiment/embedcells embedded
%   on; METHOD ("tsne", "umap", "phate", "metaviz") describes one
%   embedding, from what SC_TSNE, SC_UMAP and SC_PHATE do with the inputs
%   EMBEDCELLS gave them.
%
%   NAME-VALUE ARGUMENTS
%     NumHVGs        genes kept by the ranker; 0 = all genes were used
%     NumGenes       genes in the matrix embedded
%     Ranker         "analytic", "splinefit" or "brennecke" (PKG.I_SELECTHVGS)
%     NumAppended    marker genes appended to the HVGs
%     MarkerWeight   their weight, 1 = not up-weighted
%     NumDims        embedding dimensions
%     NumCells       cells embedded
%
%   See also SINGLECELLEXPERIMENT/EMBEDCELLS, PKG.I_METHODSREPORT.

arguments
    part (1,1) string
    options.NumHVGs (1,1) double = 0
    options.NumGenes (1,1) double = NaN
    options.Ranker (1,1) string = "analytic"
    options.NumAppended (1,1) double = 0
    options.MarkerWeight (1,1) double = 1
    options.NumDims (1,1) double = 2
    options.NumCells (1,1) double = NaN
end

refkeys = "cai2020";
weighted = options.MarkerWeight > 1 && options.NumAppended > 0;
switch part
    case "features"
        if options.NumHVGs > 0
            switch options.Ranker
                case "analytic"
                    how = "ranked on raw counts by their distance from the " + ...
                        "mean–variability–dropout curve expected under a " + ...
                        "gamma-Poisson model, a closed-form version of the " + ...
                        "spline-DV curve (Gatlin et al., 2025)";
                    refkeys = [refkeys, "gatlin2025"];
                case "splinefit"
                    how = "ranked on raw counts by their distance from a " + ...
                        "spline fitted to mean, coefficient of variation and " + ...
                        "dropout rate (spline-DV; Gatlin et al., 2025)";
                    refkeys = [refkeys, "gatlin2025"];
                otherwise
                    how = "ranked on raw counts by their excess coefficient " + ...
                        "of variation over technical noise (Brennecke et " + ...
                        "al., 2013)";
                    refkeys = [refkeys, "brennecke2013"];
            end
            paragraph = "Cells were embedded on the " + pkg.i_thousands(options.NumHVGs) + ...
                " most variable genes, " + how + ".";
        else
            paragraph = "Cells were embedded on all genes, without " + ...
                "feature selection.";
        end
        if options.NumAppended > 0
            paragraph = paragraph + " " + pkg.i_thousands(options.NumAppended) + ...
                " marker genes from PanglaoDB were added to them.";
            refkeys = [refkeys, "franzen2019"];
        end

    case "tsne"
        prep = i_prep(weighted, options.MarkerWeight);
        if options.NumGenes > 500
            npc = 50;
            if options.NumDims >= 10, npc = 5*options.NumDims; end
            paragraph = sprintf("A %d-dimensional t-SNE embedding (van der " + ...
                "Maaten and Hinton, 2008) was computed with MATLAB's tsne " + ...
                "(Barnes-Hut, default perplexity of 30) on the first %d " + ...
                "principal components of ", options.NumDims, npc) + prep + ".";
        else
            paragraph = sprintf("A %d-dimensional t-SNE embedding (van der " + ...
                "Maaten and Hinton, 2008) was computed with MATLAB's tsne " + ...
                "(exact algorithm, default perplexity of 30) on ", ...
                options.NumDims) + prep + ".";
        end
        refkeys = [refkeys, "vandermaaten2008"];

    case "umap"
        prep = i_prep(weighted, options.MarkerWeight);
        if options.NumCells > 500
            input = "the first 50 principal components (randomized PCA) of " + prep;
        else
            input = prep;
        end
        if isMATLABReleaseOlderThan('R2026a')
            impl = "the UMAP implementation of MATLAB File Exchange " + ...
                "submission 71902 (15 neighbours)";
        else
            impl = "MATLAB's umap (15 neighbours, Euclidean distance, 200 " + ...
                "epochs, reproducible mode)";
        end
        paragraph = sprintf("A %d-dimensional UMAP embedding (McInnes et " + ...
            "al., 2018) was computed with ", options.NumDims) + impl + ...
            " on " + input + ".";
        refkeys = [refkeys, "mcinnes2018"];

    case "phate"
        paragraph = sprintf("A %d-dimensional PHATE embedding (Moon et al., " + ...
            "2019) was computed with the authors' MATLAB implementation " + ...
            "(k = 10, t = 20) on up to 100 principal components of " + ...
            "library-size-normalized, square-root-transformed expression.", ...
            options.NumDims);
        refkeys = [refkeys, "moon2019"];

    otherwise
        paragraph = sprintf("A %d-dimensional %s embedding was computed.", ...
            options.NumDims, upper(part));
end
end


function prep = i_prep(weighted, weight)
% What SC_TSNE and SC_UMAP embed: EMBEDCELLS normalises first itself only
% when it up-weights markers, and then tells them not to again.
prep = "library-size-normalized (10,000 counts per cell, over the " + ...
    "selected genes), log(x + 1)-transformed expression";
if weighted
    prep = prep + sprintf(", with the added marker genes multiplied by %g", weight);
end
end
