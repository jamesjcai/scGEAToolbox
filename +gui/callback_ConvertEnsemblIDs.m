function [requirerefresh] = callback_ConvertEnsemblIDs(src)
%CALLBACK_CONVERTENSEMBLIDS Replace Ensembl gene IDs in SCE.G with gene symbols.
%
%   Some datasets -- h5ad files from public archives especially -- carry
%   Ensembl stable IDs in SCE.G rather than gene symbols, which makes every
%   marker-based tool in the GUI useless until they are converted.
%
%   Only names that look like Ensembl gene IDs are touched, so running this
%   on a mixed gene list, or twice in a row, is safe. An ID with no entry in
%   the mapping table keeps its original name rather than being blanked.

requirerefresh = false;

[FigureHandle, sce] = gui.gui_getfigsce(src);

g = string(sce.g);
[isEns, bare] = i_findensemblids(g);
nEns = sum(isEns);

if nEns == 0
    gui.myHelpdlg(FigureHandle, ...
        ['None of the gene names look like Ensembl IDs, so there is ', ...
        'nothing to convert.'], '');
    return;
end

species = i_guessspecies(bare(isEns), FigureHandle);
if species == "", return; end   % user cancelled the species prompt

answer = gui.myQuestdlg(FigureHandle, ...
    sprintf(['%d of %d gene names look like Ensembl IDs. Convert them ', ...
    'to %s gene symbols?'], nEns, numel(g), species), ...
    'Convert Ensembl IDs');
if ~strcmp(answer, 'Yes'), return; end

% The undo snapshot is taken by the menu handler in scgeatoolApp.mlapp,
% before this runs, which is the convention for every other destructive
% operation and what TESTS/UNDOGUARDTEST asserts. SingleCellExperiment is a
% handle class, so it has to happen before anything here writes to SCE.G.
idx = find(isEns);
symbols = string(pkg.e_ensembl2symbol(bare(idx), char(species)));

% E_ENSEMBL2SYMBOL returns its input untouched for an ID it cannot map, so
% an unchanged entry means "not found". Those keep the name the dataset came
% with, version suffix and all, rather than the stripped-down ID: the point
% is to add information, never to quietly damage a name that could not be
% improved.
mapped = symbols ~= bare(idx);
newg = g;
newg(idx(mapped)) = symbols(mapped);

nMapped = sum(mapped);
if nMapped == 0
    gui.myHelpdlg(FigureHandle, ...
        sprintf(['None of the %d Ensembl IDs were found in the %s ', ...
        'mapping table. Gene names are unchanged; check that the ', ...
        'species is right.'], nEns, species), '');
    return;
end

sce.g = newg;
gui.myGuidata(FigureHandle, sce, src);
requirerefresh = true;

gui.myHelpdlg(FigureHandle, i_summary(nMapped, nEns, newg), '');
end


function [isEns, bare] = i_findensemblids(g)
%I_FINDENSEMBLIDS Flag Ensembl gene IDs and strip any version suffix.
%   Stable IDs are a species prefix, a G for "gene", then digits --
%   ENSG00000121410 for human, ENSMUSG00000064336 for mouse. Many files
%   append a version, ENSG00000121410.5, which is not in the mapping table
%   and has to come off before lookup.
%
%   Transcript and protein IDs (ENST..., ENSP...) are deliberately not
%   matched: they are not what the gene-level table is keyed on.

g = string(g);
isEns = matches(g, regexpPattern("ENS[A-Z]*G\d+(\.\d+)?"));
bare = extractBefore(g, caseInsensitivePattern(".") | textBoundary("end"));
bare(ismissing(bare)) = g(ismissing(bare));
end


function species = i_guessspecies(ids, FigureHandle)
%I_GUESSSPECIES Infer species from the ID prefix, asking only when unsure.
%   ENSG is human and ENSMUSG is mouse; those two cover the mapping tables
%   shipped with the toolbox. Anything else, or a mix, has to be asked.

nHuman = sum(startsWith(ids, "ENSG"));
nMouse = sum(startsWith(ids, "ENSMUSG"));

if nHuman > 0 && nMouse == 0
    species = "human";
    return;
end
if nMouse > 0 && nHuman == 0
    species = "mouse";
    return;
end

answer = gui.myQuestdlg(FigureHandle, ...
    'Which species mapping table should be used?', ...
    'Convert Ensembl IDs', {'human', 'mouse'}, 'human');
if isempty(answer) || ~ismember(answer, {'human', 'mouse'})
    species = "";
    return;
end
species = string(answer);
end


function msg = i_summary(nMapped, nEns, newg)
%I_SUMMARY Describe what changed, including any newly duplicated symbols.
%   Distinct Ensembl IDs can share a symbol, so a conversion that succeeds
%   can still leave two rows with the same name. That is worth saying: gene
%   lookups elsewhere in the GUI take the first match.

msg = sprintf("Converted %d of %d Ensembl IDs to gene symbols.", ...
    nMapped, nEns);

nUnmapped = nEns - nMapped;
if nUnmapped > 0
    msg = msg + sprintf(newline + ...
        "%d ID(s) were not in the mapping table and keep their original name.", ...
        nUnmapped);
end

nDup = numel(newg) - numel(unique(newg));
if nDup > 0
    msg = msg + sprintf(newline + ...
        "%d gene name(s) are now duplicated, because separate Ensembl IDs " + ...
        "share a symbol.", nDup);
end
end
