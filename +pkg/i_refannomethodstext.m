function txt = i_refannomethodstext(ref, options)
%I_REFANNOMETHODSTEXT Methods paragraph for reference-based cell type annotation.
%
%   pkg.i_refannomethodstext(ref, Name=Value)        prints it and, given SCE, keeps it
%   txt = pkg.i_refannomethodstext(ref, Name=Value)  returns the paragraph instead
%
%   Writes what SC_CELLTYPEANNOREF did with reference REF, at its default
%   settings, in a form a user can paste into a manuscript's Methods. The
%   Allen Brain Cell Atlas reference (SC_ALLENBRAINREF) is described with
%   how it was built and cited; any other reference by its ref.Source and
%   number of types. Which tests are described follows the rule
%   SC_CELLTYPEANNOREF uses: a test runs only if Method asks for it and REF
%   has what it needs (Profile for the correlation, Markers for the
%   overlap test). Printing goes through PKG.I_METHODSREPORT, which adds
%   the references and, when SCE is given, stores the paragraph as the
%   "annotation" step for SC_METHODSTEXT. Called by SC_ANNOTATECELLS.
%
%   NAME-VALUE ARGUMENTS
%     Method        "both" (default), "correlation" or "markers", as passed
%                   to SC_CELLTYPEANNOREF
%     NumClusters   clusters annotated; NaN leaves the number out
%     NumGenesUsed  genes the correlation used; NaN leaves it out
%     NumAgree      clusters on which both tests named the same type
%                   ("both"), or that got a type at all ("markers");
%                   NaN leaves it out
%     SCE           the SingleCellExperiment annotated, to keep the text
%                   with it; [] (default) prints only
%     Print         false to store in SCE without printing (default true)
%
%   See also SC_CELLTYPEANNOREF, SC_ALLENBRAINREF, PKG.I_CTAMETHODSTEXT.

arguments
    ref (1,1) struct
    options.Method (1,1) string {mustBeMember(options.Method, ...
        ["both", "correlation", "markers"])} = "both"
    options.NumClusters (1,1) double = NaN
    options.NumGenesUsed (1,1) double = NaN
    options.NumAgree (1,1) double = NaN
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

hasProfile = isfield(ref, "Profile") && ~isempty(ref.Profile);
hasMarkers = isfield(ref, "Markers") && ~isempty(ref.Markers);
useCorrelation = hasProfile && options.Method ~= "markers";
useMarkers = hasMarkers && options.Method ~= "correlation";

[refText, buildText, isAllen] = in_describeref(ref);
refkeys = ["acri2026", "cai2020"];
if isAllen
    refkeys = [refkeys, "yao2023"];
end

paragraph = "Cell types were assigned to the " + clusters + ...
    " by comparison with " + refText + ", following SAHA (Acri et al., " + ...
    "2026) as implemented in " + toolbox + " (Cai, 2020). The query was " + ...
    "not integrated with the reference: each cluster was compared with " + ...
    "per-type summary statistics only." + buildText;

if useCorrelation
    if isnan(options.NumGenesUsed)
        genes = "";
    else
        genes = "; " + pkg.i_thousands(options.NumGenesUsed) + " used here";
    end
    paragraph = paragraph + " Counts were scaled to library size and " + ...
        "log(x + 1) transformed, and each cluster's mean profile was " + ...
        "compared with each reference type's mean profile over the shared " + ...
        "genes whose means differed between the query clusters (one-way " + ...
        "ANOVA, Benjamini–Hochberg-adjusted P < 0.01; at most the 2,000 " + ...
        "with the largest F" + genes + "). Each gene was z-scored across " + ...
        "the groups of its own dataset, so that differences in depth and " + ...
        "platform cancel, and clusters and types were compared by Spearman " + ...
        "correlation. Each cluster was assigned its best-correlated type; " + ...
        "types within 0.05 of the best were re-decided, as in SingleR " + ...
        "(Aran et al., 2019), by Euclidean distance over the genes whose " + ...
        "z-scores separated them by at least one.";
    refkeys = [refkeys, "aran2019"];
end

if useMarkers
    if useCorrelation
        lead = " Independently, each";
    else
        lead = " Each";
    end
    paragraph = paragraph + lead + " cluster's positive markers (Wilcoxon " + ...
        "rank-sum test against all other cells, Bonferroni-adjusted " + ...
        "P < 0.05, fold change > 1.5, detected in ≥ 50% of the cluster's " + ...
        "cells and < 25% of other cells; at most 100) were tested for " + ...
        "over-representation among each type's markers by a hypergeometric " + ...
        "test over the shared genes, with Benjamini–Hochberg correction " + ...
        "over all cluster–type pairs.";
    refkeys = [refkeys, "benjamini1995"];
end

count = "";
if ~isnan(options.NumAgree) && ~isnan(options.NumClusters)
    count = sprintf(" (%d of %d clusters)", options.NumAgree, ...
        options.NumClusters);
end
if useCorrelation && useMarkers
    paragraph = paragraph + " The label was taken from the correlation, " + ...
        "and was considered supported where the marker test named the " + ...
        "same type at FDR < 0.05" + count + ".";
elseif useMarkers
    paragraph = paragraph + " Each cluster was assigned the type with " + ...
        "the smallest FDR if it was below 0.05" + count + ", and was " + ...
        "labelled ""undetermined"" otherwise.";
end

if nargout > 0
    txt = paragraph;
else
    pkg.i_methodsreport(options.SCE, "annotation", paragraph, refkeys, ...
        Print=options.Print);
end
end


function [text, buildText, isAllen] = in_describeref(ref)
% The Allen reference is recognised by its citation, which SC_ALLENBRAINREF
% sets and a hand-made reference would not.
numTypes = numel(ref.Types);
source = "";
buildText = "";
if isfield(ref, "Source")
    source = string(ref.Source);
end
isAllen = isfield(ref, "Citation") && ...
    contains(string(ref.Citation), "s41586-023-06812-z");
if isAllen
    if endsWith(source, "class") && ~endsWith(source, "subclass")
        level = "class";
        levels = "classes";
    else
        level = "subclass";
        levels = "subclasses";
    end
    text = "the Allen Brain Cell Atlas whole mouse brain taxonomy " + ...
        "(CCN20230722; Yao et al., 2023) at the " + level + " level (" + ...
        numTypes + " " + levels + ")";
    buildText = " Each " + level + "'s reference profile was the " + ...
        "cell-weighted mean of the published per-cluster mean expression " + ...
        "(log2(CPM + 1)), over the genes some " + level + " expressed at " + ...
        "≥ 1, and its markers were the genes whose mean was ≥ 1 and " + ...
        "exceeded the cell-weighted mean of all other " + levels + " by " + ...
        "≥ 1 (two-fold), at most 100.";
    return
end
if strlength(source) > 0
    text = "a reference of " + numTypes + " cell types (" + source + ")";
else
    text = "a user-supplied reference of " + numTypes + " cell types";
end
end
