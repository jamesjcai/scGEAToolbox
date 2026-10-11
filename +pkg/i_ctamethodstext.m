function txt = i_ctamethodstext(options)
%I_CTAMETHODSTEXT Methods paragraph for a CTA-scoring cell type annotation.
%
%   pkg.i_ctamethodstext(Name=Value)        prints it and, given SCE, keeps it
%   txt = pkg.i_ctamethodstext(Name=Value)  returns the paragraph instead
%
%   Writes what RUN.ML_ALONA did, from the settings it ran with, in a form a
%   user can paste into a manuscript's Methods. Printing goes through
%   PKG.I_METHODSREPORT, which adds the references and, when SCE is given,
%   stores the paragraph as the "annotation" step for SC_METHODSTEXT.
%   Printed once by each entry point - SC_CELLTYPEANNO, the GUI's database
%   marker annotation, SingleCellExperiment/assigncelltype and the CLI - and
%   never by RUN.ML_ALONA itself, which the GUI calls once per cluster.
%
%   NAME-VALUE ARGUMENTS
%     Species      "human" or "mouse": which PanglaoDB marker sheet
%     Summary      "median" (default, as the paper) or "mean"
%     Alpha        FDR cut-off for the Unknown rule (0.05)
%     NumClusters  clusters annotated; NaN leaves the number out
%     NumUnknown   clusters labelled Unknown; NaN leaves it out
%     Manual       true when the user picked each cluster's type from the
%                  ranked candidates, so the Unknown rule did not apply
%     SCE          the SingleCellExperiment annotated, to keep the text
%                  with it; [] (default) prints only
%     Print        false to store in SCE without printing (default true)
%
%   See also RUN.ML_ALONA, SC_CELLTYPEANNO.

arguments
    options.Species (1,1) string = "human"
    options.Summary (1,1) string {mustBeMember(options.Summary, ["median", "mean"])} = "median"
    options.Alpha (1,1) double = 0.05
    options.NumClusters (1,1) double = NaN
    options.NumUnknown (1,1) double = NaN
    options.Manual (1,1) logical = false
    options.SCE = []
    options.Print (1,1) logical = true
end

version = string(pkg.i_get_versionnum());
if strlength(version) > 0
    toolbox = "scGEAToolbox " + version;
else
    toolbox = "scGEAToolbox";
end
if isnan(options.NumClusters)
    clusters = "clusters";
else
    clusters = sprintf("%d clusters", options.NumClusters);
end

if options.Summary == "median"
    summary = "summarised per cluster by each gene's median expression";
else
    summary = "summarised per cluster by each gene's mean expression " + ...
        "(the original method uses the median, which was kept for calling " + ...
        "a gene expressed)";
end

if options.Manual
    assignment = "For each cluster, the cell type was then chosen by hand " + ...
        "from the ten highest-scoring candidates.";
else
    assignment = "Each cluster was assigned its highest-scoring cell type " + ...
        "if that type's markers were over-represented among the genes " + ...
        "expressed in the cluster (median > 0; one-sided Fisher's exact " + ...
        "test, Benjamini–Hochberg-adjusted P ≤ " + string(options.Alpha) + ...
        "), and was labelled ""Unknown"" otherwise";
    if ~isnan(options.NumUnknown) && ~isnan(options.NumClusters)
        assignment = assignment + sprintf(" (%d of %d clusters)", ...
            options.NumUnknown, options.NumClusters);
    end
    assignment = assignment + ".";
end

paragraph = "Cell types were assigned to the " + clusters + ...
    " by cell-type activity (CTA) scoring (Franzén et al., 2019) " + ...
    "against the " + lower(options.Species) + " marker sets of " + ...
    "PanglaoDB, as implemented in " + toolbox + " (Cai, 2020). Counts " + ...
    "were scaled to library size, log2(x + 1) transformed and " + summary + ...
    ". Each marker gene was weighted by its specificity, " + ...
    "w = 1 + √((max f − f) / (max f − min f)), where f is the number of " + ...
    "cell types it marks, and the CTA score of a cell type was the " + ...
    "weighted sum of its markers' expression divided by the cube root of " + ...
    "its number of markers. " + assignment;

if nargout > 0
    txt = paragraph;
else
    pkg.i_methodsreport(options.SCE, "annotation", paragraph, ...
        ["franzen2019", "cai2020"], Print=options.Print);
end
end
