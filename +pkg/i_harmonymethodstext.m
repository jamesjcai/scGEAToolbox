function [paragraph, refkeys] = i_harmonymethodstext(options)
%I_HARMONYMETHODSTEXT Methods text for Harmony batch correction.
%
%   [paragraph, refkeys] = pkg.i_harmonymethodstext(Name=Value)
%
%   Describes what GUI.CALLBACK_HARMONY (native, RUN.ML_HARMONY),
%   GUI.CALLBACK_HARMONYR or GUI.CALLBACK_HARMONYPY did.
%
%   NAME-VALUE ARGUMENTS
%     Backend     "native" (default), "r" or "python"
%     UsePCs      true: corrected the PCs of PKG.E_CELLPCS and re-embedded
%                 with UMAP; false: corrected the embedding on screen
%     NumBatches  batches corrected across
%     NumDims     dimensions of the embedding corrected or produced
%
%   See also RUN.ML_HARMONY, PKG.E_CELLPCS, PKG.I_METHODSREPORT.

arguments
    options.Backend (1,1) string {mustBeMember(options.Backend, ["native", "r", "python"])} = "native"
    options.UsePCs (1,1) logical = true
    options.NumBatches (1,1) double = NaN
    options.NumDims (1,1) double = 2
end

refkeys = ["korsunsky2019", "cai2020"];
if isnan(options.NumBatches)
    across = "across batches";
else
    across = sprintf("across the %d batches", options.NumBatches);
end
switch options.Backend
    case "native"
        impl = "the native MATLAB implementation in scGEAToolbox (theta = 2, " + ...
            "ridge lambda = 1, up to 10 iterations, min(100, max(2, N/30)) " + ...
            "soft clusters for N cells, with the kernel bandwidth set by " + ...
            "perplexity calibration rather than the published fixed 0.1)";
    case "r"
        impl = "the harmony R package with its default parameters";
    otherwise
        impl = "harmonypy with its default parameters";
end
paragraph = "Batch effects " + across + " were corrected with Harmony " + ...
    "(Korsunsky et al., 2019), using " + impl;
if options.UsePCs && options.Backend == "native"
    paragraph = paragraph + ", on the top 30 principal components of the " + ...
        "2,000 most variable genes (library-size-normalized, " + ...
        "log(x + 1)-transformed and scaled per gene, clipped at 10). The " + ...
        sprintf("corrected components were embedded in %d dimensions with ", ...
        options.NumDims) + "UMAP (McInnes et al., 2018) and used for " + ...
        "subsequent graph-based clustering.";
    refkeys = [refkeys, "mcinnes2018"];
else
    paragraph = paragraph + sprintf(", applied directly to the " + ...
        "%d-dimensional embedding.", options.NumDims);
end
end
