function [tmpcelltypev, speciestag] = callback_Brush4Celltypes(src)
%CALLBACK_BRUSH4CELLTYPES Assign a cell type to brushed cells (one-time analysis).
%   [tmpcelltypev, speciestag] = callback_Brush4Celltypes(src) determines the
%   cell type of the currently brushed cells, then labels the representative
%   cell with a datatip. Labels are not written to the sce, which
%   GUI.I_CONFIRMONCE says once a session rather than on every selection.
%
%   The markers are PanglaoDB's by default. The other choice is a customized
%   marker list - typed, loaded from a file, or started from ScTypeDB.
%   GUI.I_GETCUSTOMMARKERS remembers it for the session, so the next selection
%   opens on the markers just used and takes one OK to confirm them.
%
%   The brush is offered a chance to grow to the whole group it sits in first,
%   through GUI.I_EXPANDBRUSHED - the same question Delete Selected Cells asks.
%
%   src can be the scgeatool app or a component of a plain figure.
%
%   Outputs are returned so an app caller can persist them across calls:
%       [app.tmpcelltypev, app.speciestag] = gui.callback_Brush4Celltypes(app);
%   Both outputs are pre-populated from src, so an early return leaves the
%   caller's existing state unchanged.
%
%   See also gui.i_getcustommarkers, gui.i_sessionmarkers,
%   pkg.i_celltypebrushed, pkg.i_markerweights, pkg.e_determinecelltype.

[FigureHandle, sce] = gui.gui_getfigsce(src);

isapp = isa(src, 'matlab.apps.AppBase');
if isapp
    speciestag = src.speciestag;
    tmpcelltypev = src.tmpcelltypev;
    h = src.h;
else
    speciestag = [];
    tmpcelltypev = [];
    h = findobj(FigureHandle, 'Type', 'Scatter');
    if ~isscalar(h) && ~isempty(h)
        h = h(1);
    end
end

if ~pkg.i_isvalid(h)
    gui.myWarndlg(FigureHandle, 'No plot available.');
    return;
end
ptsSelected = logical(h.BrushData.');
if ~any(ptsSelected)
    gui.myHelpdlg(FigureHandle, ...
        "No cells are selected. Please use the data brush tool to " + ...
        "select cells for cell type assignment.");
    return;
end
% Once per session. The notice is worth reading the first time and is in the
% way of the twentieth selection, which is how this tool is used.
if ~gui.i_confirmonce(FigureHandle, 'brush4celltypes:notsaved', ...
        ['This is a one-time analysis. Cell type labels will not be ' ...
        'saved to the data. You will not be asked again this session.'], ...
        'Assign Cell Type to Selected Cells')
    return;
end

% A few cells of a cluster are how a cluster gets brushed, and the type of
% those few is rarely the question being asked. GUI.I_EXPANDBRUSHED offers
% the whole group instead, the way Delete Selected Cells does, and asks only
% when there is something to ask: the brush has to sit inside one group of
% several. Everything after this works on what it returns, including the
% cell the datatip ends up on.
[ptsSelected, letdoit] = gui.i_expandbrushed(ptsSelected, sce, FigureHandle);
if ~letdoit, return; end

% Which markers to score the selection against. The customized path needs no
% species: the genes are the user's own, so there is nothing to look up.
answer = gui.myQuestdlg(FigureHandle, ['Which marker genes should the ' ...
    'selected cells be scored against?'], 'Marker Genes', ...
    {'PanglaoDB (built-in)', 'Customized marker genes...'}, ...
    'PanglaoDB (built-in)');
if isempty(answer), return; end
usecustom = strcmp(answer, 'Customized marker genes...');

if usecustom
    Tm = gui.i_getcustommarkers(FigureHandle, sce);
    if isempty(Tm), return; end
    [wvalu, wgene, celltypev, markergenev] = pkg.i_markerweights(Tm);
else
    if isempty(speciestag)
        speciestag = gui.i_selectspecies(2, false, FigureHandle, speciestag);
    end
    if isempty(speciestag), return; end
end

fw = gui.myWaitbar(FigureHandle);
if usecustom
    % The same scoring GUI.CALLBACK_DETERMINECELLTYPECLUSTERS and
    % SC_CSUBTYPEANNO use for a user's own markers.
    [Tct] = pkg.e_determinecelltype(sce, ptsSelected, wvalu, wgene, ...
        celltypev, markergenev);
else
    [Tct] = pkg.i_celltypebrushed(sce.X, sce.g, sce.s, ptsSelected, ...
        speciestag, "all", "panglaodb", false);
end
ctxt = Tct.C1_Cell_Type;
gui.myWaitbar(FigureHandle, fw);

if gui.i_isuifig(FigureHandle)
    [indx, tf] = gui.myListdlg(FigureHandle, ctxt, 'Select cell type');
else
    [indx, tf] = listdlg('PromptString', ...
        {'Select cell type'}, 'SelectionMode', 'single', ...
        'ListString', ctxt, 'ListSize', [220, 300]);
end
if tf ~= 1, return; end
ctxt = Tct.C1_Cell_Type{indx};

ctxt = strrep(ctxt, '_', '\_');
delete(findall(FigureHandle, 'Type', 'hggroup'));
if isempty(tmpcelltypev) || length(tmpcelltypev) ~= sce.NumCells
    tmpcelltypev = cell(sce.NumCells, 1);
end

siv = sce.s(ptsSelected, :);
si = mean(siv, 1);
[k] = dsearchn(siv, si);
idx = find(ptsSelected);
if isempty(k) || k > length(idx)
    gui.myWarndlg(FigureHandle, 'Could not find nearest cell.');
    return;
end
tmpcelltypev{idx(k)} = ctxt;

% h is a handle, so updating the template here updates the caller's plot.
if pkg.i_isvalid(h)
    rowx = dataTipTextRow('', tmpcelltypev);
    h.DataTipTemplate.DataTipRows = rowx;
    datatip(h, 'DataIndex', idx(k));
end
end
