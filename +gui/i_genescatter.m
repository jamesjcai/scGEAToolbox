function i_genescatter(T, parentfig)

if nargin < 2, parentfig = []; end
hx = gui.myFigure(parentfig);

g = T.gene;

hFig = hx.FigHandle;
hAx = hx.AxHandle;
if isempty(hAx), hAx = gca; end

hx.addCustomButton('off', @in_callback_HighlightGenes, 'plotpicker-qqplot.gif', 'Highlight top DV genes');
hx.addCustomButton('off', @in_callback_HighlightSelectedGenes, 'curve-array.jpg', 'Highlight selected genes');
hx.addCustomButton('off', @in_callback_ExportGeneNames, 'bookmark-book.jpg', 'Export Selected HVG gene names...');
hx.addCustomButton('off', @in_callback_ExportTable, 'floppy-disk-arrow-in.jpg', 'Export HVG Table...');
hx.addCustomButton('off', @in_callback_EnrichrHVGs, 'plotpicker-andrewsplot.gif', 'Enrichment analysis...');
hx.addCustomButton('off', @in_callback_ChangeAlphaValue, 'Brightness-3--Streamline-Core.jpg', 'Change MarkerFaceAlpha value');

isHighlighting = false;
h = scatter(hAx, T.de_coef, T.dv_coef, 'filled');
xlabel(hAx, 'DE Coefficient');
ylabel(hAx, 'DV Coefficient');

if ~isempty(g)
    dt = datacursormode(hFig);
    dt.UpdateFcn = {@in_myupdatefcn3, g};
end
xline(hAx, 0,'r');
yline(hAx, 0,'r');
hx.show;

function in_callback_ChangeAlphaValue(~, ~)
if h.MarkerFaceAlpha <= 0.05
    h.MarkerFaceAlpha = 1;
else
    h.MarkerFaceAlpha = h.MarkerFaceAlpha - 0.1;
end
end

function in_callback_HighlightSelectedGenes(~,~)
[glist] = gui.i_selectngenes(string(T.gene), [], hFig);
if ~isempty(glist)
    [y,idx]=ismember(glist, T.gene);
    idx=idx(y);
    % idv = zeros(1, length(hvgidx));
    % idv(idx)=1;
    % h.BrushData = idv;
    for k=1:length(idx)
        dt = datatip(h,'DataIndex',idx(k));
    end
end
end

function in_callback_HighlightGenes(~, ~)

% The toolbar button stays live while the dialog is up, so without this
% guard a second click re-enters the callback and opens another dialog.
if isHighlighting, return; end
isHighlighting = true;
restoreFlag = onCleanup(@() in_clearhighlightflag());

[~, hvgidx] = sort(T.dv_pval);

idx = zeros(1, height(T));
h.BrushData = idx;
% Parent the dialog on this figure, not on the caller's window: a dialog
% modal to parentfig leaves the figure owning the button clickable.
k = gui.i_inputnumk(200, 1, 2000, [], hFig);
if isempty(k), return; end
k = min(k, numel(hvgidx));
idx(hvgidx(1:k)) = 1;
h.BrushData = idx;
end

function in_clearhighlightflag()
isHighlighting = false;
end

function in_callback_ExportTable(~, ~)
gui.i_exporttable(T, true, 'Tmementores', 'MementoRsTable', [], [], hFig);
% Tdegenelist
% 'Tviolindata','ViolinPlotTable'
% 'Thvgreslist', 'HVGResultTable'
end

function in_callback_ExportGeneNames(~, ~)
ptsSelected = logical(h.BrushData.');
if ~any(ptsSelected)
    gui.myWarndlg(hFig, "No gene is selected.");
    return;
end
fprintf('%s selected.\n', pkg.i_plural(sum(ptsSelected), 'gene'));

gselected=g(ptsSelected);
[yes,idx]=ismember(gselected, g);
Tx=T(idx,:);
Tx=sortrows(Tx,'dv_pval','ascend');
if ~all(yes), error('Running time error.'); end
tgenes=Tx.gene;


labels = {'Save selected gene names to variable:',...
    'Save HVG table:'};
vars = {'g','T'};
values = {tgenes,T};
export2wsdlg(labels, vars, values, ...
    'Save Data to Workspace');
end

function in_callback_EnrichrHVGs(~, ~)
ptsSelected = logical(h.BrushData.');
if ~any(ptsSelected)
    gui.myWarndlg(hFig,"No gene is selected.");
    return;
end
fprintf('%s selected.\n', pkg.i_plural(sum(ptsSelected), 'gene'));

gselected=g(ptsSelected);
[yes,idx]=ismember(gselected, T.gene);
Tx=T(idx,:);
Tx=sortrows(Tx,'dv_pval','ascend');
if ~all(yes), error('Running time error.'); end
tgenes=Tx.gene;

gui.i_enrichtest(tgenes, g, numel(tgenes));
end


function txt = in_myupdatefcn3(src, event_obj, g)

if isequal(get(src, 'Parent'), hAx)
    idx = event_obj.DataIndex;
    txt = g(idx);
else
    txt = num2str(event_obj.Position(2));
end
end

end
