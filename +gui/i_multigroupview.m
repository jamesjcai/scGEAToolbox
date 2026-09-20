function hFig = i_multigroupview(sce, cv, panelTitles, parentfig, figname)
%I_MULTIGROUPVIEW One embedding panel per grouping, with the panels linked.
%
%   hFig = gui.i_multigroupview(sce, cv, panelTitles)
%   hFig = gui.i_multigroupview(sce, cv, panelTitles, parentfig, figname)
%
%   CV is a cell array of per-cell grouping vectors, each with sce.NumCells
%   elements - cluster ids, cell types, a stashed annotation, anything
%   GUI.I_SELECTNSTATES returns. PANELTITLES names each one. Every panel draws
%   the SAME embedding (sce.s) coloured by a different grouping, which is what
%   makes the three links between them meaningful:
%
%     rotate    one camera for all panels, so a 3-D view stays comparable
%     brush     cells selected in one panel light up in every other, which is
%               how "which cells do these two disagree about" gets answered
%     datatip   the label under the cursor, per panel
%
%   plus a toolbar button that writes each group's name at its centroid.
%
%   Returns the figure handle, or [] if the cells have no embedding to plot
%   on.
%
%   This is the drawing half of GUI.CALLBACK_MULTIGROUPINGVIEW, lifted out so
%   that a caller holding its own set of groupings - annotation history, say -
%   gets the same figure without the picker in front of it.
%
%   See also GUI.CALLBACK_MULTIGROUPINGVIEW, GUI.I_SELECTNSTATES,
%   GUI.CALLBACK_COMPARECELLTYPEANNOTATIONS, GUI.I_GSCATTER3.

if nargin < 5, figname = ''; end
if nargin < 4, parentfig = []; end

hFig = [];
if isempty(cv), return; end
if ~iscell(cv), cv = {cv}; end
if nargin < 3 || isempty(panelTitles)
    panelTitles = "Grouping " + string(1:numel(cv));
end
panelTitles = string(panelTitles);

% PKG.E_HASEMBEDDING first, and the shape tests after it. The shape tests
% alone could not fire: the SingleCellExperiment constructor fills S with
% randn(nCells, 3) when the caller supplies none, so size(sce.s, 1) ==
% sce.NumCells and size(sce.s, 2) == 3 hold on an object that has never been
% embedded, and the panels then get drawn over random Gaussian coordinates.
if ~pkg.e_hasembedding(sce) || size(sce.s, 1) ~= sce.NumCells || ...
        size(sce.s, 2) < 2
    gui.myWarndlg(parentfig, ['The cells have no 2-D embedding to plot on. ', ...
        'Run an embedding first.']);
    return;
end

npanel = numel(cv);
for k = 1:npanel
    if numel(cv{k}) ~= sce.NumCells
        error('gui:i_multigroupview:sizeMismatch', ...
            ['Grouping %d has %d elements and the dataset has %d cells. ', ...
            'Every grouping must give one value per cell.'], ...
            k, numel(cv{k}), sce.NumCells);
    end
end

hx = gui.myFigure(parentfig);
hFig = hx.FigHandle;
if strlength(figname) > 0
    hFig.Name = figname;
    hFig.NumberTitle = 'off';
end
hFig.Position(3) = hFig.Position(3)*1.8;

% An explicit layout rather than bare NEXTTILE, which falls back to whichever
% figure is current and grows a flow layout there.
tl = tiledlayout(hFig, 'flow', 'TileSpacing', 'compact', 'Padding', 'compact');

axesv = cell(npanel, 1);
hv = cell(npanel, 1);
cmapv = cell(npanel, 1);
labelv = cell(npanel, 1);
for k = 1:npanel
    labelv{k} = string(cv{k}(:));
    axesv{k} = nexttile(tl);
    hv{k} = gui.i_gscatter3(sce.s, cv{k}, 1, 1, axesv{k});
    title(axesv{k}, panelTitles(k), 'Interpreter', 'none');
    cmapv{k} = colormap(axesv{k});
end

hx.addCustomButton('off', @in_showgrouplabels, "label.jpg", ...
    "Show group labels");
hx.show(parentfig);

dt = datacursormode(hFig);
dt.UpdateFcn = @in_datatip;

% The link object has to outlive this function or the cameras come apart
% again the moment it is garbage collected. It used to be built inside an
% EVALIN('base', ...) whose result was assigned to nothing at all, so the
% panels have not in fact rotated together; keeping it on the figure fixes
% that and takes the base workspace out of it.
setappdata(hFig, 'MultiGroupCameraLink', ...
    linkprop(vertcat(axesv{:}), {'CameraPosition', 'CameraUpVector'}));

rotate3d(hFig, 'on');
hBr = brush(hFig);
hBr.ActionPostCallback = @in_onbrush;

% COLORMAP on a tile can reset its neighbours, so the per-panel maps are put
% back after every panel exists.
for k = 1:npanel
    colormap(axesv{k}, cmapv{k});
end


    function in_onbrush(~, event)
        % Mirror the brushed cells into every other panel. Same cells, same
        % positions, so the selection means the same thing in each.
        srcidx = in_findaxes(event.Axes);
        if isempty(srcidx), return; end
        d = hv{srcidx}.BrushData;
        for kb = 1:npanel
            if kb ~= srcidx
                hv{kb}.BrushData = d;
            end
        end
    end


    function txt = in_datatip(target, event)
        txt = '';
        tipidx = in_findaxes(ancestor(target, 'axes'));
        if isempty(tipidx), return; end
        lbl = labelv{tipidx};
        txt = gui.i_escapeunderscore(lbl(event.DataIndex));
    end


    function in_showgrouplabels(~, ~)
        % Toggle: labels on if none are showing, otherwise clear them.
        hastip = false;
        for kl = 1:npanel
            tips = findobj(hv{kl}, 'Type', 'datatip');
            if ~isempty(tips)
                delete(tips);
                hastip = true;
            end
        end
        if hastip, return; end

        for kl = 1:npanel
            [gidx, gname] = findgroups(labelv{kl});
            if max(gidx) >= in_maxlabelgroups(), continue; end
            gname = gui.i_escapeunderscore(gname);
            hv{kl}.DataTipTemplate.DataTipRows = dataTipTextRow('', gname(gidx));
            for kg = 1:max(gidx)
                idx = find(gidx == kg);
                siv = sce.s(idx, :);
                % The cell nearest the group's centroid, by direct distance.
                % DSEARCHN answers the same question for a single query
                % point, but triangulates first, which on a dataset of this
                % size is where the button appeared to hang.
                [~, near] = min(sum((siv - mean(siv, 1)).^2, 2));
                datatip(hv{kl}, 'DataIndex', idx(near));
            end
        end
    end


    function idx = in_findaxes(ax)
        idx = [];
        if isempty(ax), return; end
        for ka = 1:npanel
            if isequal(ax, axesv{ka})
                idx = ka;
                return;
            end
        end
    end

end


function n = in_maxlabelgroups()
% Past this many groups the centroid labels overlap into an unreadable mass,
% and the panel is better read by brushing and hovering.
n = 50;
end
