function [paragraph, refkeys] = i_dvmethodstext(options)
%I_DVMETHODSTEXT Methods text for a differential variability run.
%
%   [paragraph, refkeys] = pkg.i_dvmethodstext(Name=Value)
%
%   Describes what SC_DVG did. With no permutations its p-values are, by
%   SC_DVG's own account, a ranking aid rather than a calibrated test, and
%   the text says so rather than present them as one.
%
%   NAME-VALUE ARGUMENTS
%     Method          "splinefit" (default), "analytic" or "brennecke"
%     Direction       "mean" (default) or "deviation"
%     NumPermutations 0, or the permutations used for the p-values
%     Group1, Group2  group names; N1, N2  cells after filtering
%     QCFiltered      each group was QC-filtered on its own first
%                     (GUI.CALLBACK_DVGENE2GROUPS)
%     NumGenesRemoved uninformative genes removed first; NaN = not stated
%
%   See also SC_DVG, PKG.I_METHODSREPORT.

arguments
    options.Method (1,1) string {mustBeMember(options.Method, ["splinefit", "analytic", "brennecke"])} = "splinefit"
    options.Direction (1,1) string {mustBeMember(options.Direction, ["mean", "deviation"])} = "mean"
    options.NumPermutations (1,1) double = 0
    options.Group1 (1,1) string = "group 1"
    options.Group2 (1,1) string = "group 2"
    options.N1 (1,1) double = NaN
    options.N2 (1,1) double = NaN
    options.QCFiltered (1,1) logical = false
    options.NumGenesRemoved (1,1) double = NaN
end

refkeys = "cai2020";
g1 = """" + options.Group1 + """";
g2 = """" + options.Group2 + """";
if ~isnan(options.N1), g1 = g1 + " (" + pkg.i_thousands(options.N1) + " cells)"; end
if ~isnan(options.N2), g2 = g2 + " (" + pkg.i_thousands(options.N2) + " cells)"; end

prep = "";
if ~isnan(options.NumGenesRemoved) && options.NumGenesRemoved > 0
    prep = prep + sprintf("%d mitochondrial, ribosomal, hemoglobin or " + ...
        "unapproved-symbol genes were removed first. ", options.NumGenesRemoved);
end
if options.QCFiltered
    prep = prep + "Each group was quality-filtered on its own (at least " + ...
        "1,000 counts and 200 detected genes per cell, under 15% " + ...
        "mitochondrial counts, genes detected in at least 15 cells). ";
end

switch options.Method
    case "splinefit"
        how = "with spline-DV (Gatlin et al., 2025): in each group, every " + ...
            "gene was placed by its log mean, log coefficient of variation " + ...
            "and dropout rate after library-size normalization, a spline " + ...
            "curve was fitted through the genes, and a gene's DV score was " + ...
            "the distance between its two deviation vectors from the curves " + ...
            "of the two groups";
        refkeys = [refkeys, "gatlin2025"];
    case "analytic"
        how = "with the spline-DV framework (Gatlin et al., 2025), using in " + ...
            "place of the fitted spline the closed-form curve of mean, " + ...
            "coefficient of variation and dropout rate expected under a " + ...
            "gamma-Poisson model; a gene's DV score was the distance between " + ...
            "its two deviation vectors from the curves of the two groups";
        refkeys = [refkeys, "gatlin2025"];
    otherwise
        how = "as the difference between the two groups in each gene's " + ...
            "excess squared coefficient of variation over the technical-" + ...
            "noise fit of Brennecke et al. (2013)";
        refkeys = [refkeys, "brennecke2013"];
end
paragraph = prep + "Differential variability between " + g1 + " and " + ...
    g2 + " was assessed " + how + ".";

if options.Direction == "mean"
    paragraph = paragraph + " Genes were labelled up or down by whether " + ...
        "their library-size-normalized mean was higher in the first group.";
else
    paragraph = paragraph + " Genes were labelled up or down by whether " + ...
        "they deviated further from the expected curve in the first group.";
end
if options.NumPermutations > 0
    paragraph = paragraph + sprintf(" P values came from %d random " + ...
        "reassignments of the pooled cells to two groups of the original " + ...
        "sizes, each gene's distance scaled by the root mean square of its " + ...
        "own permuted distances and compared with the pooled scaled " + ...
        "permuted distances.", options.NumPermutations);
else
    paragraph = paragraph + " Genes were ranked by DV score; the " + ...
        "accompanying P values, from a half-normal null scaled by the " + ...
        "median score over all genes, were used for ranking only, as they " + ...
        "are not calibrated.";
end
end
