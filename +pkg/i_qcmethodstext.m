function paragraph = i_qcmethodstext(options)
%I_QCMETHODSTEXT Methods paragraph for cell and gene quality control.
%
%   paragraph = pkg.i_qcmethodstext(Name=Value)
%
%   Describes what SC_QCFILTER did (through SingleCellExperiment's qcfilter
%   and qcfilterwhitelist) or, with Iterative=false and LibSize=NaN, what
%   CLI.CMD_FILTER did, from the thresholds actually used. Values below 1
%   are described as SC_QCFILTER reads them: MinCells as a fraction of
%   cells, LibSize as a quantile of library size.
%
%   NAME-VALUE ARGUMENTS
%     LibSize        minimum total counts per cell; NaN = not applied
%     MtRatio        mitochondrial fraction cut-off
%     MtInclusive    true: cells AT or above MtRatio go (SC_QCFILTER);
%                    false: only cells above it (CMD_FILTER)
%     MinCells       a gene must be detected in at least this many cells
%     MinGenes       a cell must have at least this many detected genes
%     Iterative      filters repeated until nothing more is removed
%     NumWhitelist   genes kept regardless of the filters
%     RemovedMtGenes mitochondrial genes were dropped afterwards
%     CellsBefore, CellsAfter, GenesBefore, GenesAfter   the counts
%
%   See also SC_QCFILTER, PKG.I_METHODSREPORT.

arguments
    options.LibSize (1,1) double = 1000
    options.MtRatio (1,1) double = 0.15
    options.MtInclusive (1,1) logical = true
    options.MinCells (1,1) double = 15
    options.MinGenes (1,1) double = 200
    options.Iterative (1,1) logical = true
    options.NumWhitelist (1,1) double = 0
    options.RemovedMtGenes (1,1) logical = false
    options.CellsBefore (1,1) double = NaN
    options.CellsAfter (1,1) double = NaN
    options.GenesBefore (1,1) double = NaN
    options.GenesAfter (1,1) double = NaN
end

pct = @(x) string(round(100*x, 2)) + "%";
num = @pkg.i_thousands;

if options.MtInclusive
    mt = pct(options.MtRatio) + " or more";
else
    mt = "more than " + pct(options.MtRatio);
end
first = "Cells in which mitochondrial genes (symbols beginning with " + ...
    """mt-"", in either case) accounted for " + mt + " of the counts " + ...
    "were removed first.";

if options.MinCells < 1
    genes = "genes detected in fewer than " + pct(options.MinCells) + ...
        " of cells";
else
    genes = "genes detected in fewer than " + num(options.MinCells) + " cells";
end
cellcrit = "fewer than " + num(options.MinGenes) + " detected genes";
if ~isnan(options.LibSize)
    if options.LibSize < 1
        lib = "total counts below the " + string(round(100*options.LibSize, 2)) + ...
            "th percentile of library size";
    else
        lib = "fewer than " + num(options.LibSize) + " total counts";
    end
    cellcrit = lib + " or " + cellcrit;
end
if options.Iterative
    second = "Then, repeating until no further cells or genes were " + ...
        "removed, " + genes + " were removed, and cells with " + ...
        cellcrit + " were removed.";
else
    second = "Then " + genes + " were removed, and cells with " + ...
        cellcrit + " were removed; no library-size threshold was applied.";
    if ~isnan(options.LibSize)
        second = "Then " + genes + " were removed, and cells with " + ...
            cellcrit + " were removed.";
    end
end

paragraph = "Low-quality cells and rarely detected genes were removed " + ...
    "from the raw count matrix with scGEAToolbox (Cai, 2020). " + ...
    first + " " + second;
if options.NumWhitelist > 0
    paragraph = paragraph + " " + num(options.NumWhitelist) + ...
        " genes of interest were kept regardless of these criteria.";
end
if options.RemovedMtGenes
    paragraph = paragraph + " Mitochondrial genes were then removed.";
end
if ~any(isnan([options.CellsBefore, options.CellsAfter, ...
        options.GenesBefore, options.GenesAfter]))
    paragraph = paragraph + " This retained " + num(options.CellsAfter) + ...
        " of " + num(options.CellsBefore) + " cells and " + ...
        num(options.GenesAfter) + " of " + num(options.GenesBefore) + " genes.";
end
paragraph = paragraph + " Counts were not normalized at this step.";
end

