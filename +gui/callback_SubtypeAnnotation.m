function [requirerefresh] = callback_SubtypeAnnotation(src, usecustom)
% CALLBACK_SUBTYPEANNOTATION - Annotate subtypes of an annotated cell type.
%
%   requirerefresh = gui.callback_SubtypeAnnotation(src)
%   requirerefresh = gui.callback_SubtypeAnnotation(src, true)
%
% With the bundled subtype markers (the default), this offers the primary cell
% types that both the data and assets/PanglaoDB/cellsubtypes.xlsx know about,
% then hands each to SC_CSUBTYPEANNO, which isolates those cells, re-embeds and
% re-clusters them on their own, scores the new clusters against the subtype
% markers, and merges the subtype labels back.
%
% Which labels count as a primary type is decided by PKG.I_MATCHPRIMARYTYPE, so
% the menu is available for a gross annotation spelled "T cells_{3}",
% "CD8+ T cells" or "T memory cells", not only for a literal "T cells".
%
% With USECUSTOM true the marker list comes from the user instead - typed,
% loaded from a file, or started from the bundled subtypes and edited. That
% lifts the two limits of the bundled path: the population to subdivide is
% picked from the labels actually in the data rather than from the 20-odd
% primary types the table covers, and its subtypes are whatever the marker list
% names. Everything after that is the same run.
%
% see also: sc_csubtypeanno, gui.i_getcustomsubtypemarkers,
%           pkg.i_matchprimarytype, pkg.i_subtypeoverlap,
%           pkg.i_subtypecandidates, gui.i_updateannotatemenu,
%           gui.callback_CollapseCellSubtypes

requirerefresh = false;

% The second argument used to be the menu event, which is ignored. Only a
% logical scalar means "use customized markers"; anything else is an event.
if nargin < 2 || ~(islogical(usecustom) && isscalar(usecustom))
    usecustom = false;
end

if isa(src, "SingleCellExperiment")
    sce = src;
    FigureHandle = [];
else
    [FigureHandle, sce] = gui.gui_getfigsce(src);
end

% "undetermined" is what the constructor fills in, not an annotation.
if ~pkg.i_hascelltypelabels(sce)
    gui.myErrordlg(FigureHandle, ['Cells are not annotated yet. Annotate ' ...
        'cell types first, then come back for their subtypes.']);
    return;
end
labels = string(sce.c_cell_type_tx);

if usecustom
    [celltypetarget_list, selections, Tsub] = ...
        in_pickcustom(FigureHandle, sce, labels);
else
    [celltypetarget_list, selections, Tsub] = ...
        in_pickbundled(FigureHandle, sce, labels);
end
if isempty(celltypetarget_list), return; end

answer = gui.myQuestdlg(FigureHandle, 'How to label cell type with subtype', ...
    'Choose format', {'Type (Subtype)','Type_{Subtype}','Subtype'}, ...
    'Type (Subtype)');
if isempty(answer), return; end
switch answer
    case 'Type (Subtype)'
        formatid = 1;
    case 'Type_{Subtype}'
        formatid = 2;
    case 'Subtype'
        formatid = 0;
end

% Re-embedding and re-clustering a subset takes a while, and doing it for
% several types in a row takes several times as long.
% The species decides which primary marker file SC_CSUBTYPEANNO appends
% from, and its SPECIESTAG argument used to be ignored. Guessing from the
% symbol casing is what PKG.I_GUESSSPECIES is for; it had no callers.
speciestag = pkg.i_guessspecies(sce.g);

% Keep the labels this run is about to overwrite, the way the other
% annotation handlers do. It matters more here than there: with the 'Subtype'
% format the new label keeps no trace of the parent type, so without this the
% main cell type is simply gone and GUI.CALLBACK_COLLAPSECELLSUBTYPES has nothing
% to strip.
stashname = pkg.i_stashcelltypehistory(sce);

fw = gui.myWaitbar(FigureHandle);
a.FigureHandle = FigureHandle;
a.fw = fw;
for k = 1:length(celltypetarget_list)
    fw = gui.myWaitbar(FigureHandle, fw, false, '', ...
        sprintf('Annotating subtypes of %s...', celltypetarget_list(k)), ...
        (k-1)/length(celltypetarget_list));
    opts = struct();
    if usecustom
        opts.SubtypeMarkers = Tsub;
        opts.CellSelection = selections{k};
    end
    [sce] = sc_csubtypeanno(sce, celltypetarget_list(k), formatid, speciestag, a, opts);
end
gui.myWaitbar(FigureHandle, fw);
gui.myGuidata(FigureHandle, sce, src);
requirerefresh = true;

msg = sprintf('%s now in the data.', pkg.i_plural( ...
    numel(unique(string(sce.c_cell_type_tx))), 'cell type label'));
if formatid ~= 0
    msg = msg + " Annotate > Collapse Cell Subtypes back to Cell Type collapses them back.";
end
gui.myHelpdlg(FigureHandle, msg + gui.i_stashnotice(stashname));
end

function [ctypelist, selections, Tsub] = in_pickbundled(FigureHandle, sce, labels)
% The types cellsubtypes.xlsx has markers for, as before. SELECTIONS stays
% empty: SC_CSUBTYPEANNO finds the cells itself through PKG.I_MATCHPRIMARYTYPE.

ctypelist = strings(0, 1);
selections = {};
Tsub = [];

% The same question the Annotate menu asks before enabling this item.
[candidates, matched, primarytypes, resolved] = pkg.i_subtypecandidates(sce);
if isempty(candidates)
    gui.myErrordlg(FigureHandle, sprintf(['None of the %d cell type label(s) ' ...
        'in your data is a cell type that cellsubtypes.xlsx has subtype ' ...
        'markers for. Supported: %s. Use "Annotate Cell Subtypes Using ' ...
        'Customized Marker Genes..." to subdivide any other cell type.'], ...
        numel(unique(labels)), strjoin(primarytypes, ', ')));
    return;
end

% Say how many cells each choice covers, and under how many different labels:
% the count is what tells the user that "T cells" here also takes in the cells
% annotated as "CD8+ T cells". Cells that already carry a subtype are counted
% apart, because SC_CSUBTYPEANNO leaves them where they are and a single total
% would promise to re-annotate cells it will not touch.
items = strings(size(candidates));
for k = 1:numel(candidates)
    isk = matched == candidates(k);
    items(k) = sprintf('%s (%s, %s)', candidates(k), ...
        pkg.i_plural(sum(isk & ~resolved), 'cell'), ...
        pkg.i_plural(numel(unique(labels(isk & ~resolved))), 'label'));
    if any(isk & resolved)
        items(k) = items(k) + sprintf(' - %d more already subtyped', ...
            sum(isk & resolved));
    end
end

prompt = ['Cells of the selected type are isolated, re-embedded and ' ...
    're-clustered on their own before their subtypes are annotated.'];
if gui.i_isuifig(FigureHandle)
    [indx2, tf2] = gui.myListdlg(FigureHandle, cellstr(items), ...
        'Subtype Annotation', [], true, true, [420, 260], prompt);
else
    [indx2, tf2] = listdlg('PromptString', {prompt}, ...
        'SelectionMode', 'multiple', 'ListString', ...
        cellstr(items), 'ListSize', [420, 260]);
end

if tf2 ~= 1, return; end
ctypelist = candidates(indx2);
end

function [ctypelist, selections, Tsub] = in_pickcustom(FigureHandle, sce, labels)
% Any label in the data, and a marker list the user supplies. Several labels
% can be picked at once and are then subdivided together as one population:
% "T cells_{1}" and "T cells_{4}" are usually two clusters of the same thing,
% and splitting them apart before re-clustering would only put the same
% boundary back.

ctypelist = strings(0, 1);
selections = {};
Tsub = [];

[ulabels, ~, back] = unique(labels(:));
counts = accumarray(back, 1);
[counts, ord] = sort(counts, 'descend');
ulabels = ulabels(ord);

items = arrayfun(@(k) sprintf('%s (%s)', ulabels(k), ...
    pkg.i_plural(counts(k), 'cell')), (1:numel(ulabels))', ...
    'UniformOutput', false);

prompt = ['Select the cell type(s) to subdivide. Those cells are pooled, ' ...
    'isolated, re-embedded and re-clustered on their own, and each new ' ...
    'cluster is labelled by its best match among the subtype markers you ' ...
    'supply next.'];
if gui.i_isuifig(FigureHandle)
    [indx, tf] = gui.myListdlg(FigureHandle, items, ...
        'Subtype Annotation', [], true, true, [420, 300], prompt);
else
    [indx, tf] = listdlg('PromptString', {prompt}, ...
        'SelectionMode', 'multiple', 'ListString', items, ...
        'ListSize', [420, 300]);
end
if tf ~= 1 || isempty(indx), return; end

picked = ulabels(indx);
targetname = in_targetname(picked);

Tsub = gui.i_getcustomsubtypemarkers(FigureHandle, sce, targetname);
if isempty(Tsub), return; end

gui.i_warnmissingmarkers(FigureHandle, sce, Tsub, ...
    'Customized Subtype Markers');

ctypelist = targetname;
selections = {ismember(labels, picked)};
end

function [name] = in_targetname(picked)
% What to call the pooled population in "Type (Subtype)". The cluster index
% SCE.ASSIGNCELLTYPE appends is dropped, so picking "T cells_{1}" and
% "T cells_{4}" gives "T cells" rather than a two-part name.

base = unique(erase(picked, "_{" + digitsPattern + "}"), 'stable');
if isscalar(base)
    name = base;
else
    name = strjoin(base, "+");
end
end
