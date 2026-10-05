function [phase, position, T] = r_tricycle(X, genelist, options)
%R_TRICYCLE  Cell cycle position and phase with the R package tricycle.
%   [phase, position] = run.r_tricycle(X, genelist) projects each cell onto
%   tricycle's cell-cycle space, learned from mouse neurosphere data, and
%   returns its position on the cycle, an angle in [0, 2*pi), with a G1, S
%   or G2M label read off it.
%
%   [phase, position, T] = run.r_tricycle(...) also returns a table with
%   the position and the two projection coordinates (pc1, pc2) per cell.
%
%   Name-value arguments:
%   Species        : "human" or "mouse". Guessed from the gene-symbol case
%                    when omitted.
%   ReferenceX     : raw count matrix (genes x cells) of a reference sample,
%                    such as untreated control cells.
%   ReferenceGenes : gene names for the rows of ReferenceX.
%   WorkDir        : folder for the exchange files. Defaults to a temporary
%                    folder.
%
%   tricycle centres every gene on the mean of the cells it scores before
%   projecting them, so, like Seurat's module scores, a sample dominated by
%   one phase moves its own origin and every angle with it. With ReferenceX
%   the sample and the reference are normalised together and the angle is
%   taken around the reference cells' mean instead. Only genes present in
%   both are used, matched exactly and then ignoring case.
%
%   Labels: S from 0.5*pi, G2M from pi to 1.25*pi, G1 from 1.25*pi round to
%   0.5*pi. The S and G2M starts are the tricycle vignette's landmarks; the
%   vignette runs G2M on to 1.75*pi, but on the example data cells past
%   1.25*pi have no G2M module-score signal, so they are counted as
%   M/early G1 and labelled G1. The labels are approximate; the position is
%   the primary output.
%
%   Zheng et al. 2022, Genome Biology. Bioconductor package "tricycle".
%
%   See also SC_CELLCYCLESCORE, RUN.R_SEURATCELLCYCLE

arguments
    X
    genelist
    options.Species {mustBeTextScalar} = ""
    options.ReferenceX = []
    options.ReferenceGenes = []
    options.WorkDir = pkg.i_tempdirfile()
end

genelist = string(genelist(:));
if size(X, 1) ~= numel(genelist)
    error("run:r_tricycle:size", ...
        "X has %d rows but GENELIST has %d names. Pass one gene name per row of X.", ...
        size(X, 1), numel(genelist));
end
species = lower(string(options.Species));
if strlength(species) == 0
    species = string(pkg.i_guessspecies(genelist));
end
if ~ismember(species, ["human", "mouse"])
    error("run:r_tricycle:species", ...
        "tricycle has references for human and mouse only, not ""%s"".", species);
end

useReference = ~isempty(options.ReferenceX);
if useReference
    [X, genelist, refX] = i_aligntoreference(X, genelist, ...
        options.ReferenceX, options.ReferenceGenes);
end

oldpth = pwd();
cleanupCwd = onCleanup(@() cd(oldpth));
[isok, msg, codepath] = commoncheck_R('R_tricycle');
if ~isok, error('%s', msg); end
if ~isempty(options.WorkDir) && isfolder(options.WorkDir), cd(options.WorkDir); end

tmpfilelist = {'input.h5', 'reference.h5', 'species.txt', 'output.csv'};
% Always cleared first: a stale reference.h5 would turn reference mode on,
% and a stale output.csv would be read back as this run's result.
pkg.i_deletefiles(tmpfilelist);
cleanupFiles = onCleanup(@() pkg.i_deletefiles(tmpfilelist));
pkg.e_writeh5(sparse(double(X)), genelist, 'input.h5');
if useReference
    pkg.e_writeh5(sparse(double(refX)), genelist, 'reference.h5');
end
writelines(species, 'species.txt');

Rpath = getpref('scgeatoolbox', 'rexecutablepath', []);
if isempty(Rpath)
    error('R environment has not been set up.');
end
pkg.i_runrcode(fullfile(codepath, 'script.R'), Rpath);

if ~isfile('output.csv')
    error('run:r_tricycle:noOutput', ...
        ['R finished but did not write output.csv to %s. The R console ', ...
        'output above should say why.'], pwd);
end
T = readtable('output.csv', 'ReadVariableNames', true);
if height(T) ~= size(X, 2)
    error('run:r_tricycle:outputSize', ...
        'tricycle returned %d cells for %d in the input.', height(T), size(X, 2));
end
position = T.position;
phase = i_position2phase(position);
T.phase = phase;
end

function phase = i_position2phase(position)
% S and G2M start where the tricycle vignette puts them (0.5*pi, pi). G2M
% ends at 1.25*pi, not the vignette's 1.75*pi: past 1.25*pi the G2M module
% score is back at background, so those cells -- late M, early G1 -- are G1.
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
    error('run:r_tricycle:referenceSize', ...
        ['ReferenceX has %d rows but ReferenceGenes has %d names. ', ...
        'Pass one gene name per row of ReferenceX.'], size(refX, 1), numel(refGenes));
end
[found, loc] = ismember(genelist, refGenes);
[foundCI, locCI] = ismember(upper(genelist(~found)), upper(refGenes));
loc(~found) = locCI;
found(~found) = foundCI;
if ~any(found)
    error('run:r_tricycle:referenceNoGenes', ...
        'None of the genes in GENELIST are in ReferenceGenes. Check that both use the same gene symbols.');
end
X = X(found, :);
genelist = genelist(found);
refX = refX(loc(found), :);
end
