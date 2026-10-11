function [phase, position, T] = ml_tricycle(X, genelist, options)
%ML_TRICYCLE  Cell cycle position and phase with tricycle, in MATLAB.
%   [phase, position] = run.ml_tricycle(X, genelist) projects each cell onto
%   tricycle's cell-cycle space, learned from mouse neurosphere data, and
%   returns its position on the cycle, an angle in [0, 2*pi), with a G1, S
%   or G2M label read off it.
%
%   [phase, position, T] = run.ml_tricycle(...) also returns a table with
%   the position, the two projection coordinates (pc1, pc2) and the phase
%   per cell.
%
%   A MATLAB port of run.r_tricycle, with the same inputs and outputs and
%   no R. It uses tricycle's own reference matrix, copied unchanged into
%   assets/Tricycle (GPL-3, see the README there), and reproduces
%   project_cycle_space and estimate_cycle_position:
%     1. log2(counts/sizeFactor + 1), size factors = library size over its
%        mean (scuttle's logNormCounts);
%     2. keep the reference genes, matched by exact symbol, first
%        occurrence -- mouse symbols, or tricycle's human mapping;
%     3. centre each gene and project onto the two rotation columns;
%     4. position = atan2(pc2, pc1), wrapped to [0, 2*pi).
%   Reference genes missing from X are dropped, not imputed, as tricycle
%   does.
%
%   Name-value arguments:
%   Species        : "human" or "mouse". Guessed from the gene-symbol case
%                    when omitted.
%   ReferenceX     : raw count matrix (genes x cells) of a reference sample,
%                    such as untreated control cells.
%   ReferenceGenes : gene names for the rows of ReferenceX.
%
%   Without ReferenceX, each gene is centred on the mean of the cells being
%   scored, so a sample dominated by one phase moves its own origin and
%   every angle with it. With ReferenceX the sample and the reference are
%   normalised together and each gene is centred on the reference cells'
%   mean instead. Only genes present in both are used, matched exactly and
%   then ignoring case.
%
%   Labels: S from 0.5*pi, G2M from pi to 1.25*pi, G1 from 1.25*pi round to
%   0.5*pi; see RUN.R_TRICYCLE for why. Cells with no counts get position
%   NaN and phase "undetermined".
%
%   Zheng et al. 2022, Genome Biology. Bioconductor package "tricycle".
%
%   See also RUN.R_TRICYCLE, SC_CELLCYCLECOMPARE, SC_CELLCYCLESCORE

arguments
    X {mustBeNumeric}
    genelist
    options.Species {mustBeTextScalar} = ""
    options.ReferenceX = []
    options.ReferenceGenes = []
end

genelist = string(genelist(:));
if size(X, 1) ~= numel(genelist)
    error("run:ml_tricycle:size", ...
        "X has %d rows but GENELIST has %d names. Pass one gene name per row of X.", ...
        size(X, 1), numel(genelist));
end
species = lower(string(options.Species));
if strlength(species) == 0
    species = string(pkg.i_guessspecies(genelist));
end
if ~ismember(species, ["human", "mouse"])
    error("run:ml_tricycle:species", ...
        "tricycle has references for human and mouse only, not ""%s"".", species);
end

nCells = size(X, 2);
useReference = ~isempty(options.ReferenceX);
if useReference
    [X, genelist, refX] = i_aligntoreference(X, genelist, ...
        options.ReferenceX, options.ReferenceGenes);
    % Normalised together, so sample and reference share one size-factor scale.
    X = [X, refX];
end

ref = i_loadreference();
if species == "human"
    refSymbols = ref.human_symbol;
else
    refSymbols = ref.mouse_symbol;
end
% ismember gives the first occurrence, as R's match does.
[found, loc] = ismember(refSymbols, genelist);
if ~any(found)
    error("run:ml_tricycle:noGenes", ...
        ['None of tricycle''s %d %s reference genes are in GENELIST. ', ...
        'Check the Species and that GENELIST holds gene symbols.'], ...
        numel(refSymbols), species);
end
rotation = [ref.pc1_rot(found), ref.pc2_rot(found)];

libSize = full(sum(X, 1, "double"));
sizeFactor = libSize/mean(libSize);
Y = full(double(X(loc(found), :)))./sizeFactor;
Y = log2(Y + 1);

if useReference
    centre = mean(Y(:, nCells+1:end), 2, "omitnan");
    Y = Y(:, 1:nCells);
else
    centre = mean(Y, 2, "omitnan");
end
embedding = (Y - centre)'*rotation;

position = mod(atan2(embedding(:, 2), embedding(:, 1)), 2*pi);
phase = i_position2phase(position);
T = table(position, embedding(:, 1), embedding(:, 2), phase, ...
    VariableNames=["position", "pc1", "pc2", "phase"]);
end

function ref = i_loadreference()
persistent cached
if isempty(cached)
    f = fullfile(fileparts(fileparts(mfilename("fullpath"))), ...
        "assets", "Tricycle", "neuroRef.tsv");
    opts = detectImportOptions(f, FileType="text", Delimiter="\t");
    opts = setvartype(opts, ["mouse_symbol", "human_symbol", "ensembl"], "string");
    opts = setvartype(opts, ["pc1_rot", "pc2_rot"], "double");
    cached = readtable(f, opts);
end
ref = cached;
end

function phase = i_position2phase(position)
% Same boundaries as run.r_tricycle.
phase = repmat("G1", numel(position), 1);
phase(position >= 0.5*pi & position < pi) = "S";
phase(position >= pi & position < 1.25*pi) = "G2M";
phase(isnan(position)) = "undetermined";
end

function [X, genelist, refX] = i_aligntoreference(X, genelist, refX, refGenes)
% Keep the genes in both, in GENELIST order: exact names first, so two
% genes that differ only in case keep their own rows, then ignoring case.
refGenes = string(refGenes(:));
if size(refX, 1) ~= numel(refGenes)
    error("run:ml_tricycle:referenceSize", ...
        ['ReferenceX has %d rows but ReferenceGenes has %d names. ', ...
        'Pass one gene name per row of ReferenceX.'], size(refX, 1), numel(refGenes));
end
[found, loc] = ismember(genelist, refGenes);
[foundCI, locCI] = ismember(upper(genelist(~found)), upper(refGenes));
loc(~found) = locCI;
found(~found) = foundCI;
if ~any(found)
    error("run:ml_tricycle:referenceNoGenes", ...
        "None of the genes in GENELIST are in ReferenceGenes. Check that both use the same gene symbols.");
end
X = X(found, :);
genelist = genelist(found);
refX = refX(loc(found), :);
end
