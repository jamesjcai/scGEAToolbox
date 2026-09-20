function [T, info] = sc_fgsea(genelist, stat, opts)
%SC_FGSEA  Preranked gene set enrichment against Enrichr libraries.
%
%   T = SC_FGSEA(genelist) runs preranked GSEA on an already-ordered gene
%   list, testing it against the Enrichr libraries the fGSEA workflow has
%   always used, and returns the sets enriched toward the top of the list.
%
%   T = SC_FGSEA(genelist, stat) ranks by STAT instead of by the order of
%   GENELIST.
%
%   This is the native replacement for RUN.R_FGSEA. The enrichment score and
%   its permutation null come from SC_GSETTEST(Method="gsea"); this function
%   supplies the gene sets, the ranking convention and the filtering that
%   the R script wrapped around fgsea.
%
%   It is not fgsea's algorithm. fgsea's contribution is a multilevel
%   sampling scheme that resolves p-values far below 1/NumPerm, and that
%   resolution is not a nicety. A plain permutation p-value stops at
%   1/(NumPerm+1), so every strongly enriched set ties there: on a test
%   where one pathway was placed at the very top of the ranking it scored a
%   perfect ES of 1.0 and still came back at FDR 0.098, tied with sets of
%   NES 1.5, and nothing cleared a 0.05 cutoff at all. Here the tail is
%   extrapolated with a generalized Pareto fit instead
%   (SC_GSETTEST PValueTail="gpd"), which recovers that pathway at rank 1
%   with FDR 1e-5 and puts the runner-up at 1.3e-4 against an independently
%   measured 1e-4. It is an approximation where fgsea is exact.
%
%   USAGE:
%     Tdeg = sc_deg(X1, X2, genelist);
%     T = sc_fgsea(Tdeg.gene, Tdeg.tScore);
%
%   INPUTS:
%     genelist - gene symbols. When STAT is omitted the list is taken to be
%                ranked already, most interesting first.
%     stat     - ranking statistic, one value per gene, larger meaning
%                "further toward the top". Pass [] to use the list order.
%
%   NAME-VALUE:
%     Libraries       - Enrichr libraries to test against. Default is the
%                       five the R script used, so results stay comparable.
%     RemoveRibosomal - drop ribosomal genes before ranking (default true,
%                       as in the R path). They dominate many rankings for
%                       reasons that are rarely the biology of interest.
%     Direction       - "up" (default) scores enrichment toward the top
%                       only, matching fgsea's scoreType="pos". "both"
%                       reports depletion as well.
%     MinSize/MaxSize - set size limits after intersecting with GENELIST
%                       (defaults 15 and 500, the GSEA convention).
%     NumPerm         - permutations for the null (default 5000). The null
%                       is drawn once per distinct set size, so this costs
%                       far less than one null per set would.
%     FDRCutoff       - keep sets at or below this FDR (default 0.05). Set 1
%                       to keep everything.
%     Verbose         - print progress (default true).
%
%   OUTPUTS:
%     T    - one row per enriched set: SetName, SetSize, Stat, AUC, PValue,
%            FDR, ES, NES. Sorted by FDR.
%     info - struct with the membership matrix, the gene universe and the
%            statistic actually used.
%
%   See also SC_GSETTEST, PKG.E_GETENRICHRSETS, RUN.R_FGSEA.

arguments
    genelist
    stat double = []
    opts.Libraries string = ["KEGG_2019_Human", "BioPlanet_2019", ...
        "GO_Biological_Process_2018", "GO_Molecular_Function_2018", ...
        "Reactome_2016"]
    opts.RemoveRibosomal (1,1) logical = true
    opts.Direction (1,1) string {mustBeMember(opts.Direction, ...
        ["up", "both"])} = "up"
    opts.MinSize (1,1) double {mustBePositive, mustBeInteger} = 15
    opts.MaxSize (1,1) double {mustBePositive} = 500
    opts.NumPerm (1,1) double {mustBePositive, mustBeInteger} = 5000
    opts.FDRCutoff (1,1) double {mustBeNonnegative} = 0.05
    opts.Verbose (1,1) logical = true
end

genelist = upper(string(genelist(:)));
if isempty(genelist)
    error("sc_fgsea:noGenes", "GENELIST is empty.");
end

if isempty(stat)
    % No statistic given, so the list order is the ranking. A descending
    % ramp reproduces that ordering exactly and keeps every value positive,
    % which is what the R script's fixed template of positive values did.
    stat = double(numel(genelist):-1:1).';
else
    stat = double(stat(:));
    if numel(stat) ~= numel(genelist)
        error("sc_fgsea:statSize", ...
            "STAT has %d values but GENELIST has %d genes.", ...
            numel(stat), numel(genelist));
    end
end

keep = isfinite(stat);
if opts.RemoveRibosomal
    try
        ribosomal = upper(string(pkg.i_get_ribosomalgenes));
        keep = keep & ~ismember(genelist, ribosomal);
    catch ME
        % The list is fetched over the network; losing it weakens the
        % ranking but does not invalidate it, so carry on and say so.
        warning("sc_fgsea:ribosomalUnavailable", ...
            "Could not fetch the ribosomal gene list (%s); keeping them.", ...
            ME.message);
    end
end
genelist = genelist(keep);
stat = stat(keep);
if numel(genelist) < opts.MinSize
    error("sc_fgsea:tooFewGenes", ...
        "Only %d genes remain after filtering.", numel(genelist));
end

if opts.Verbose
    fprintf("sc_fgsea: %d genes ranked; fetching gene sets\n", numel(genelist));
end
[setmatrx, setnames, setgenes] = pkg.e_getenrichrsets(opts.Libraries, ...
    Verbose=opts.Verbose);

[T, gsInfo] = sc_gsettest(stat, genelist, setmatrx, setnames, setgenes, ...
    Method="gsea", Direction=opts.Direction, MinSize=opts.MinSize, ...
    MaxSize=opts.MaxSize, NumPerm=opts.NumPerm, PValueTail="gpd", ...
    Sort=false, Verbose=opts.Verbose);

% fgsea's scoreType="pos" reports enrichment at the top of the list only,
% and the R script then kept ES > 0 at padj < 0.05.
if opts.Direction == "up"
    T = T(T.ES > 0, :);
end
if opts.FDRCutoff < 1
    T = T(T.FDR <= opts.FDRCutoff, :);
end
T = sortrows(T, ["FDR", "PValue"]);

if opts.Verbose
    fprintf("sc_fgsea: %s enriched at FDR<=%g\n", ...
        pkg.i_plural(height(T), 'set'), opts.FDRCutoff);
end

if nargout > 1
    info = gsInfo;
    info.Genelist = genelist;
    info.Stat = stat;
    info.Libraries = opts.Libraries;
end

end
