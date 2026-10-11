function [paragraph, refkeys] = i_demethodstext(options)
%I_DEMETHODSTEXT Methods text for a two-group differential expression run.
%
%   [paragraph, refkeys] = pkg.i_demethodstext(Name=Value)
%
%   Describes what SC_DEG computes and, when Cutoffs is given, the cutoffs
%   PKG.E_PROCESSDETABLE applied to call genes up or down.
%
%   NAME-VALUE ARGUMENTS
%     Method        "wilcoxon" (SC_DEG method 1, default) or "ttest" (2)
%     Group1, Group2  group names; N1, N2  their cell counts
%     Cutoffs       USEDSET from PKG.E_PROCESSDETABLE: {mindiffpct,
%                   minabsolfc, apvaluecut, sortby, cutoffmode, rawpcut};
%                   {} when no genes were called (the CLI saves the table)
%     NumUp, NumDown  genes called
%     PerCellType   true when the comparison was repeated within each cell
%                   type (GUI.CALLBACK_DEGENE2GROUPSBATCH)
%
%   See also SC_DEG, PKG.E_PROCESSDETABLE, PKG.I_METHODSREPORT.

arguments
    options.Method (1,1) string {mustBeMember(options.Method, ["wilcoxon", "ttest"])} = "wilcoxon"
    options.Group1 (1,1) string = "group 1"
    options.Group2 (1,1) string = "group 2"
    options.N1 (1,1) double = NaN
    options.N2 (1,1) double = NaN
    options.Cutoffs cell = {}
    options.NumUp (1,1) double = NaN
    options.NumDown (1,1) double = NaN
    options.PerCellType (1,1) logical = false
end

refkeys = ["benjamini1995", "cai2020"];
g1 = """" + options.Group1 + """";
g2 = """" + options.Group2 + """";
if ~isnan(options.N1), g1 = g1 + " (" + pkg.i_thousands(options.N1) + " cells)"; end
if ~isnan(options.N2), g2 = g2 + " (" + pkg.i_thousands(options.N2) + " cells)"; end
if options.Method == "wilcoxon"
    test = "the two-sided Wilcoxon rank-sum test (with tie correction)";
else
    test = "a two-sample t-test";
end
if options.PerCellType
    scope = "Within each cell type, differential expression between " + ...
        g1 + " and " + g2;
else
    scope = "Differential expression between " + g1 + " and " + g2;
end
paragraph = scope + " was tested gene by gene with " + test + " on " + ...
    "library-size-normalized (10,000 counts per cell), log(x + 1)-" + ...
    "transformed expression. Log2 fold changes were computed as " + ...
    "log2(mean + 1) of the normalized expression in the first group " + ...
    "minus that in the second, and P values were adjusted by the " + ...
    "Benjamini–Hochberg procedure (Benjamini and Hochberg, 1995).";

c = options.Cutoffs;
if numel(c) >= 3
    calls = sprintf(" Genes were called differentially expressed at " + ...
        "adjusted P ≤ %g, |log2 fold change| ≥ %g and a difference of at " + ...
        "least %g in the fraction of cells detecting them", c{3}, c{2}, c{1});
    if numel(c) >= 6 && isnumeric(c{6}) && ~isempty(c{6}) && c{6} < 1
        calls = calls + sprintf(", with unadjusted P ≤ %g", c{6});
    end
    if numel(c) >= 5 && strcmpi(string(c{5}), "Automatic")
        calls = calls + " (cutoffs relaxed automatically until at least " + ...
            "20 genes passed in one direction)";
    end
    if ~isnan(options.NumUp) && ~isnan(options.NumDown)
        calls = calls + sprintf("; %d genes were higher and %d lower in " + ...
            "the first group", options.NumUp, options.NumDown);
    end
    paragraph = paragraph + calls + ".";
end
end
