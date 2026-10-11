function [paragraph, refkeys] = i_dpmethodstext(options)
%I_DPMETHODSTEXT Methods text for a differential gene-program run.
%
%   [paragraph, refkeys] = pkg.i_dpmethodstext(Name=Value)
%
%   Describes what SC_DPG did, including the filters it applies itself.
%
%   NAME-VALUE ARGUMENTS
%     Collection    option number of PKG.E_GETGENESETS (1 MSigDB,
%                   2 DoRothEA TF targets, 3 custom, 4 glycobiology,
%                   6 GlycoEnzOnto), or NaN
%     Species       "human" or "mouse" (MSigDB)
%     NumSetsTested gene sets scored; NumSetsReported  sets passing
%     RankNorm, BgSubtract  SC_DPG's two options
%     Group1, Group2  group names; N1, N2  their cell counts
%
%   See also SC_DPG, PKG.I_METHODSREPORT.

arguments
    options.Collection (1,1) double = NaN
    options.Species (1,1) string = ""
    options.NumSetsTested (1,1) double = NaN
    options.NumSetsReported (1,1) double = NaN
    options.RankNorm (1,1) logical = true
    options.BgSubtract (1,1) logical = true
    options.Group1 (1,1) string = "group 1"
    options.Group2 (1,1) string = "group 2"
    options.N1 (1,1) double = NaN
    options.N2 (1,1) double = NaN
end

refkeys = ["benjamini1995", "cai2020"];
g1 = """" + options.Group1 + """";
g2 = """" + options.Group2 + """";
if ~isnan(options.N1), g1 = g1 + " (" + pkg.i_thousands(options.N1) + " cells)"; end
if ~isnan(options.N2), g2 = g2 + " (" + pkg.i_thousands(options.N2) + " cells)"; end

switch options.Collection
    case 1
        sets = "MSigDB gene sets";
        if strlength(options.Species) > 0
            sets = sets + " (" + lower(options.Species) + ")";
        end
    case 2
        sets = "DoRothEA transcription-factor target sets";
    case 3
        sets = "user-defined gene sets";
    case 4
        sets = "glycobiology gene sets";
    case 6
        sets = "GlycoEnzOnto pathway gene sets";
    otherwise
        sets = "gene sets";
end
if ~isnan(options.NumSetsTested)
    sets = sprintf("%d ", options.NumSetsTested) + sets;
end

transform = "Expression was library-size-normalized (10,000 counts per " + ...
    "cell) and log(x + 1)-transformed";
if options.RankNorm
    transform = transform + "; genes were then ranked within each cell, " + ...
        "as in GSVA (Hänzelmann et al., 2013)";
    refkeys = [refkeys, "hanzelmann2013"];
end
if options.BgSubtract
    transform = transform + ", and each cell's mean was subtracted";
end

paragraph = "Gene-program activity was compared between " + g1 + " and " + ...
    g2 + " for " + sets + ". " + transform + ". The activity of a set in " + ...
    "a cell was the sum of these values over the set's genes measured in " + ...
    "the data. Activities were compared between the groups with the " + ...
    "two-sided Wilcoxon rank-sum test and P values adjusted by the " + ...
    "Benjamini–Hochberg procedure (Benjamini and Hochberg, 1995). Sets " + ...
    "with at least 5 measured genes, adjusted P < 0.01 and an absolute " + ...
    "log2 ratio of mean activities of at least 1 were reported";
if ~isnan(options.NumSetsReported)
    paragraph = paragraph + sprintf(" (%d sets)", options.NumSetsReported);
end
paragraph = paragraph + ".";
end
