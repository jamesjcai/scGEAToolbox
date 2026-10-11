function [paragraph, refkeys] = i_grnmethodstext(method, options)
%I_GRNMETHODSTEXT Methods text for gene regulatory network construction.
%
%   [paragraph, refkeys] = pkg.i_grnmethodstext(method, Name=Value)
%
%   METHOD is a key of NET.GRNMETHODS ("pcrnet", "xicor", "pearson",
%   "distcorr", "mi", "mic", "ids", "genie3", "pcrnet_denoised"), as
%   SC_GRN dispatches it.
%
%   NAME-VALUE ARGUMENTS
%     NumGenes, NumCells  size of the matrix the network was built from
%     Transform     key of GUI.I_TRANSFORMITEMS applied first, "none", or
%                   "" when the caller did not say
%     Groups        two group names when two networks were built and
%                   compared; empty for one network
%     ManifoldAlignment  the two networks were compared by scTenifoldNet's
%                   manifold alignment
%
%   See also SC_GRN, NET.GRNMETHODS, PKG.I_METHODSREPORT.

arguments
    method (1,1) string
    options.NumGenes (1,1) double = NaN
    options.NumCells (1,1) double = NaN
    options.Transform (1,1) string = ""
    options.Groups (1,:) string = strings(1, 0)
    options.ManifoldAlignment (1,1) logical = false
end

refkeys = "cai2020";
switch method
    case "pcrnet"
        how = "principal-component regression, as in scTenifoldNet (Osorio " + ...
            "et al., 2020): each gene was regressed on the first 3 principal " + ...
            "components of all other genes, and the regression coefficients " + ...
            "taken as edge weights";
        refkeys = [refkeys, "osorio2020"];
    case "pcrnet_denoised"
        how = "scTenifoldNet (Osorio et al., 2020): principal-component " + ...
            "regression networks were built on 10 random subsamples of 500 " + ...
            "cells and denoised by CANDECOMP/PARAFAC tensor decomposition";
        refkeys = [refkeys, "osorio2020"];
    case "xicor"
        how = "Chatterjee's xi correlation coefficient between every pair of " + ...
            "genes (Chatterjee, 2021)";
        refkeys = [refkeys, "chatterjee2021"];
    case "pearson"
        how = "Pearson correlation between every pair of genes, keeping " + ...
            "edges with |r| ≥ 0.3";
    case "distcorr"
        how = "distance correlation between every pair of genes (Székely " + ...
            "et al., 2007)";
        refkeys = [refkeys, "szekely2007"];
    case "mi"
        how = "binned mutual information between every pair of genes";
    case "mic"
        how = "the maximal information coefficient between every pair of " + ...
            "genes (Reshef et al., 2011), on at most 2,000 cells";
        refkeys = [refkeys, "reshef2011"];
    case "ids"
        how = "InterDependence Scores between every pair of genes " + ...
            "(Radhakrishnan et al., 2025)";
        refkeys = [refkeys, "radhakrishnan2025"];
    case "genie3"
        how = "GENIE3 (Huynh-Thu et al., 2010), regressing each gene on all " + ...
            "others with random forests and taking feature importances as " + ...
            "edge weights";
        refkeys = [refkeys, "huynhthu2010"];
    otherwise
        how = "the " + method + " method of scGEAToolbox";
end

size_ = "";
if ~isnan(options.NumGenes) && ~isnan(options.NumCells)
    size_ = " of " + pkg.i_thousands(options.NumGenes) + " genes across " + ...
        pkg.i_thousands(options.NumCells) + " cells";
elseif ~isnan(options.NumGenes)
    size_ = " of " + pkg.i_thousands(options.NumGenes) + " genes";
end
transform = i_transformtext(options.Transform);

if numel(options.Groups) == 2
    paragraph = "Gene regulatory networks" + size_ + " were constructed " + ...
        "separately for """ + options.Groups(1) + """ and """ + ...
        options.Groups(2) + """ by " + how + transform + ".";
    if options.ManifoldAlignment
        paragraph = paragraph + " The two networks were compared by " + ...
            "scTenifoldNet's manifold alignment (30 dimensions); each " + ...
            "gene's distance between its two aligned positions was " + ...
            "tested against a chi-square distribution with 1 degree of " + ...
            "freedom, and P values adjusted by the Benjamini–Hochberg " + ...
            "procedure (Benjamini and Hochberg, 1995).";
        refkeys = [refkeys, "osorio2020", "benjamini1995"];
    end
else
    paragraph = "A gene regulatory network" + size_ + " was constructed " + ...
        "by " + how + transform + ".";
end
end


function t = i_transformtext(key)
switch key
    case ""
        t = "";
    case "none"
        t = ", on the counts without transformation";
    case "libsize"
        t = ", on library-size-normalized expression (10,000 counts per cell)";
    case "log1p"
        t = ", on log(x + 1)-transformed counts";
    case "libsize_log1p"
        t = ", on library-size-normalized (10,000 counts per cell), " + ...
            "log(x + 1)-transformed expression";
    case "shiftedclr"
        t = ", on shifted centred-log-ratio-transformed expression";
    case "deseq"
        t = ", on DESeq size-factor-normalized expression";
    case "pearson_residuals"
        t = ", on analytic Pearson residuals";
    case "knn_smoothing"
        t = ", on kNN-smoothed expression";
    case "freeman_tukey"
        t = ", on Freeman-Tukey-transformed expression";
    case "magic"
        t = ", on MAGIC-imputed expression";
    case {"sctransform_r", "sctransform_matlab"}
        t = ", on SCTransform residuals";
    otherwise
        t = ", on " + key + "-transformed expression";
end
end
