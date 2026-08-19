function [tmpcelltypev, speciestag] = callback_Brush4Celltypes(src)
%CALLBACK_BRUSH4CELLTYPES Assign a cell type to brushed cells (one-time analysis).
%   [tmpcelltypev, speciestag] = callback_Brush4Celltypes(src) determines the
%   cell type of the currently brushed cells using PanglaoDB, then labels the
%   representative cell with a datatip. Labels are not written to the sce.
%
%   src can be the scgeatool app or a component of a plain figure.
%
%   Outputs are returned so an app caller can persist them across calls:
%       [app.tmpcelltypev, app.speciestag] = gui.callback_Brush4Celltypes(app);
%   Both outputs are pre-populated from src, so an early return leaves the
%   caller's existing state unchanged.

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
if ~strcmp(gui.myQuestdlg(FigureHandle, ...
        ['This is a one-time analysis. Cell type labels will ' ...
        'not be saved. Continue?']), 'Yes')
    return;
end

if isempty(speciestag)
    speciestag = gui.i_selectspecies(2, false, FigureHandle, speciestag);
end
if isempty(speciestag), return; end

fw = gui.myWaitbar(FigureHandle);
[Tct] = pkg.i_celltypebrushed(sce.X, sce.g, sce.s, ptsSelected, ...
    speciestag, "all", "panglaodb", false);
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
