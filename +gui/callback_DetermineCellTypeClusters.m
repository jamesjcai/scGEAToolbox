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

[c, cL] = findgroups(string(sce.c));

if usedefaultdb
    organtag = "all";
    databasetag = "panglaodb";
    if ~gui.gui_showrefinfo('PanglaoDB [PMID:30951143]', FigureHandle), return; end
    speciestag = gui.i_selectspecies(2, false, FigureHandle, speciestag);
    if isempty(speciestag), return; end
else
    [Tm, Tw] = pkg.i_markerlist2weight(sce, FigureHandle);
    if isempty(Tm) || isempty(Tw)
        return;
    end
    wvalu = Tw.Var2;
    wgene = string(Tw.Var1);
    celltypev = string(Tm.Var1);
    markergenev = string(Tm.Var2);
end

[manuallyselect, bestonly] = gui.i_annotemanner(FigureHandle);
if isempty(manuallyselect), return; end

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

for ix = 1:max(c)
    if ~manuallyselect
        gui.myWaitbar(FigureHandle, fw, false, '', '', ix/max(c));
    end
    ptsSelected = c == ix;

    if usedefaultdb
        [Tct] = pkg.i_celltypebrushed(sce.X, sce.g, ...
            sce.s, ptsSelected, ...
            speciestag, organtag, databasetag, bestonly);
    else
        [Tct] = pkg.e_determinecelltype(sce, ptsSelected, wvalu, ...
            wgene, celltypev, markergenev);
    end

    ctxt = Tct.C1_Cell_Type;
    if manuallyselect && length(ctxt) > 1
        if gui.i_isuifig(FigureHandle)
            [indx, tf] = gui.myListdlg(FigureHandle, ctxt, 'Select cell type');
        else
            [indx, tf] = listdlg('PromptString', {'Select cell type'}, ...
                'SelectionMode', 'single', 'ListString', ctxt, 'ListSize', [220, 300]);
        end
        if tf ~= 1, return; end
        ctxt = Tct.C1_Cell_Type{indx};
    else
        ctxt = Tct.C1_Cell_Type{1};
    end

    ctxt_raw = ctxt;
    ctxt = sprintf('%s_{%d}', ctxt_raw, ix);
    cL(ix) = ctxt;

    if ~isempty(h)
        ctxtdisp = strrep(ctxt_raw, '_', '\_');
        ctxtdisp = sprintf('%s_{%d}', ctxtdisp, ix);
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
