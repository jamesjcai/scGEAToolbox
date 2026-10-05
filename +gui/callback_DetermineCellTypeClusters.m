function [requirerefresh, speciestag] = callback_DetermineCellTypeClusters(src, usedefaultdb, livedatatips)

requirerefresh = false;
speciestag = [];

if nargin < 2, usedefaultdb = true; end
if nargin < 3, livedatatips = true; end

[FigureHandle, sce] = gui.gui_getfigsce(src);

if isa(src, 'matlab.apps.AppBase')
    speciestag = src.speciestag;
end

% Warn before replacing existing labels, and do it here rather than after
% the species and marker-database dialogs have been answered. The labels
% being replaced are stashed below, so this is about the active annotation
% changing, not about losing it.
if ~gui.i_confirmoverwritecelltype(FigureHandle, sce), return; end

% Annotation labels the groups on screen. A freshly loaded object has one
% group, the constructor's all-ones, and annotating that gives every cell
% the same type. Offer to cluster first, the same way Single Click Solution
% does, rather than letting the user find out from the result. It is asked
% here but run after the remaining dialogs, so cancelling one of them leaves
% the object as it was.
clusterfirst = false;
if numel(unique(sce.c)) < 2
    answer = gui.myQuestdlg(FigureHandle, ...
        ['All cells are in one group, so they would all get the same ' ...
        'cell type. Cluster them first (Louvain on principal ' ...
        'components, resolution 0.8)?'], '', ...
        {'Cluster First', 'Annotate As Is', 'Cancel'}, 'Cluster First');
    switch answer
        case 'Cluster First'
            clusterfirst = true;
        case 'Annotate As Is'
            % One group, one label: what the user asked for.
        otherwise
            return;
    end
end

if usedefaultdb
    organtag = "all";
    databasetag = "panglaodb";
    if ~gui.gui_showrefinfo('PanglaoDB [PMID:30951143]', FigureHandle), return; end
    speciestag = gui.i_selectspecies(2, false, FigureHandle, speciestag);
    if isempty(speciestag), return; end
else
    % The same getter the brush tool uses, so a list used in one is what the
    % other opens on. It also warns about markers this dataset does not
    % have, which the editor this replaces did not.
    Tm = gui.i_getcustommarkers(FigureHandle, sce);
    if isempty(Tm), return; end
    [wvalu, wgene, celltypev, markergenev] = pkg.i_markerweights(Tm);
end

[manuallyselect, bestonly] = gui.i_annotemanner(FigureHandle);
if isempty(manuallyselect), return; end

% What clustering changes, kept so a cancel further down can put it back.
% SCE is the app's handle and is clustered in place, so returning used to
% leave the new clusters on the object with the screen never redrawn.
before = struct('c', sce.c, 'c_cluster_id', sce.c_cluster_id, ...
    'struct_cell_clusterings', sce.struct_cell_clusterings, ...
    'struct_cell_reductions', sce.struct_cell_reductions);
if clusterfirst
    fw = gui.myWaitbar(FigureHandle);
    try
        sce = sce.clustercells([], 'louvainpc', true);
    catch ME
        gui.myWaitbar(FigureHandle, fw, true);
        gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
        return;
    end
    gui.myWaitbar(FigureHandle, fw);
    sce.c = sce.c_cluster_id;
end
% Numeric ids are grouped as numbers, so groups are visited 1, 2, ..., 10
% rather than "1", "10", "11", and the per-type subscripts follow that order.
if isnumeric(sce.c)
    [c, cL] = findgroups(sce.c);
    cL = string(cL);
else
    [c, cL] = findgroups(string(sce.c));
end

% Custom markers are scored against the other groups, all at once, by
% PKG.E_SCORECUSTOMMARKERS. This used to score each group on its own with
% PKG.E_DETERMINECELLTYPE, which sums median raw counts, so an abundant
% transcript present in every cell decided everything: on the mouse
% pancreas example ambient insulin made it label 6 of 10 true cell types
% "Beta cells". With a single group there is nothing to score against, so
% "Annotate As Is" keeps the old scoring.
scoreRelative = ~usedefaultdb && max(c) > 1;
if scoreRelative
    [customScore, customTypes] = pkg.e_scorecustommarkers(sce.X, sce.g, c, Tm);
end

% Set up live datatip handle if requested and available
h = [];
cLdisp = cL;
if livedatatips && isa(src, 'matlab.apps.AppBase') && ~isempty(src.h) && pkg.i_isvalid(src.h)
    h = src.h;
    dtp = findobj(h, 'Type', 'datatip');
    delete(dtp);
    hold(src.UIAxes, 'on');
end

if ~manuallyselect, fw = gui.myWaitbar(FigureHandle); end

rawTypes = strings(max(c), 1);
for ix = 1:max(c)
    if ~manuallyselect
        gui.myWaitbar(FigureHandle, fw, false, '', '', (ix - 1)/max(c));
    end
    ptsSelected = c == ix;

    if usedefaultdb
        [Tct] = pkg.i_celltypebrushed(sce.X, sce.g, ...
            sce.s, ptsSelected, ...
            speciestag, organtag, databasetag, bestonly);
    elseif scoreRelative
        Tct = in_rankedtypes(customScore(:, ix), customTypes);
    else
        [Tct] = pkg.e_determinecelltype(sce, ptsSelected, wvalu, ...
            wgene, celltypev, markergenev);
    end

    ctxt = Tct.C1_Cell_Type;
    if manuallyselect && length(ctxt) > 1
        if gui.i_isuifig(FigureHandle)
            [indx, tf] = gui.myListdlg(FigureHandle, ctxt, 'Select cell type', [], false);
        else
            [indx, tf] = listdlg('PromptString', {'Select cell type'}, ...
                'SelectionMode', 'single', 'ListString', ctxt, 'ListSize', [220, 300]);
        end
        if tf ~= 1
            in_cancel(sce, before, clusterfirst, h, src);
            return;
        end
        ctxt = Tct.C1_Cell_Type{indx};
    else
        ctxt = Tct.C1_Cell_Type{1};
    end

    % Number groups within each type (Fibroblasts_{1}, Fibroblasts_{2}, ...)
    % in the order they are annotated, rather than by group index, so the
    % subscripts of one type run 1, 2, 3 and the live datatips stay final.
    ctxt_raw = ctxt;
    rawTypes(ix) = string(ctxt_raw);
    typeOrdinal = nnz(rawTypes(1:ix) == rawTypes(ix));
    ctxt = sprintf('%s_{%d}', ctxt_raw, typeOrdinal);
    cL(ix) = ctxt;

    if ~isempty(h)
        ctxtdisp = strrep(ctxt_raw, '_', '\_');
        ctxtdisp = sprintf('%s_{%d}', ctxtdisp, typeOrdinal);
        cLdisp(ix) = ctxtdisp;
        row = dataTipTextRow('', cLdisp(c));
        h.DataTipTemplate.DataTipRows = row;
        if size(sce.s, 2) >= 2
            siv = sce.s(ptsSelected, :);
            si = mean(siv, 1);
            idx_pts = find(ptsSelected);
            [k] = dsearchn(siv, si);
            datatip(h, 'DataIndex', idx_pts(k));
        end
        drawnow;
    end
end

if ~isempty(h)
    hold(src.UIAxes, 'off');
end

if ~manuallyselect, gui.myWaitbar(FigureHandle, fw); end

% Keep the annotation this replaces as a new numbered cell attribute, the
% same way GUI.CALLBACK_RUNSCIMILARITY and SC_ANNOTATECELLS do.
stashname = pkg.i_stashcelltypehistory(sce);
sce.c_cell_type_tx = string(cL(c));

nx = length(unique(sce.c_cell_type_tx));
if nx > 1
    newtx = erase(sce.c_cell_type_tx, "_{" + digitsPattern + "}");
    if length(unique(newtx)) ~= nx
        if strcmp(gui.myQuestdlg(FigureHandle, 'Merge subclusters of same cell type?'), 'Yes')
            sce.c_cell_type_tx = newtx;
        end
    end
end

% Report what was assigned, and where the labels it replaced went. The
% pre-flight warning is a Yes/No gate that is easy to click through, so the
% attribute name is repeated here, once the name actually exists.
msg = sprintf('%d cell type(s) assigned to %d cells.', ...
    numel(unique(sce.c_cell_type_tx)), sce.NumCells);
gui.myHelpdlg(FigureHandle, msg + gui.i_stashnotice(stashname));

gui.myGuidata(FigureHandle, sce, src);
requirerefresh = true;
end


function Tct = in_rankedtypes(score, types)
% The table PKG.E_DETERMINECELLTYPE returns - C1_Cell_Type, C1_CTA_Score,
% best first - from one group's relative scores, so the automatic pick and
% the "choose manually" list read it unchanged. Types none of whose markers
% are measured are left out. When no list scores above zero the group
% matches none, and "Unknown" leads, as it did when no marker matched; the
% manual list still offers the rest beneath it.
isScored = isfinite(score);
[s, order] = sort(score(isScored), 'descend');
names = types(isScored);
names = cellstr(names(order));
if isempty(s) || s(1) <= 0
    names = [{'Unknown'}; names(:)];
    s = [0; s(:)];
end
Tct = table(names(:), s(:), 'VariableNames', {'C1_Cell_Type', 'C1_CTA_Score'});
end

function in_cancel(sce, before, clusterfirst, h, src)
%IN_CANCEL Undo what the run did before the user cancelled the manual pick.
%   Puts back the clustering done for 'Cluster First', and clears the
%   datatips drawn so far and the HOLD the loop set on the axes.
if clusterfirst
    sce.c = before.c;
    sce.c_cluster_id = before.c_cluster_id;
    sce.struct_cell_clusterings = before.struct_cell_clusterings;
    sce.struct_cell_reductions = before.struct_cell_reductions;
end
if ~isempty(h) && pkg.i_isvalid(h)
    delete(findobj(h, 'Type', 'datatip'));
    hold(src.UIAxes, 'off');
end
end
