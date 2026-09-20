function callback_Brush4Markers(src, ~)
%CALLBACK_BRUSH4MARKERS Find the genes that mark the brushed cells out.
%
%   gui.callback_Brush4Markers(src)
%
%   src can be the scgeatool app or a component of a plain figure.
%
%   The cell selection is resolved before anything else is asked. The
%   method question used to come first, so a user who had brushed nothing
%   -- the usual way to arrive here, by clicking the toolbar button or the
%   Annotate menu item before picking up the brush -- chose between lasso
%   and logistic regression, and only then was told there was no selection.

[FigureHandle, sce] = gui.gui_getfigsce(src);

if ~gui.i_installed('stats', FigureHandle), return; end

ptsSelected = i_getcellselection(FigureHandle, sce);
if isempty(ptsSelected), return; end

    switch gui.myQuestdlg(FigureHandle, 'Select method:','',...
            {'Lasso Regression','Logistic Regression 🐢 '}, ...
             'Lasso Regression')
        case 'Lasso Regression'
            uselasso=true;
        case 'Logistic Regression 🐢 '
            uselasso=false;
        otherwise
            return;
    end
i_Brush4TopMarkers(FigureHandle, sce, ptsSelected, uselasso);
end


function ptsSelected = i_getcellselection(FigureHandle, sce)
%I_GETCELLSELECTION Cells to contrast against the rest, or [] to abort.

ptsSelected = [];

axesh = FigureHandle.findobj('type', 'Axes');
if ~isscalar(axesh)
    gui.myWarndlg(FigureHandle, 'No plot available.');
    return;
end
h = axesh.findobj('type', 'Scatter');
if ~isscalar(h)
    gui.myWarndlg(FigureHandle, 'No plot available.');
    return;
end
assert(isequal(h, FigureHandle.findobj('type', 'Scatter')))

brushed = logical(h.BrushData.');

if ~any(brushed)
    answer = gui.myQuestdlg(FigureHandle, 'No cells are brushed/selected. You can select cells by a grouping variable. Continue?','');
    if ~strcmp(answer,'Yes'), return; end
    [ptsSelected] = gui.i_select1classcells(sce, false, FigureHandle);
    if isempty(ptsSelected), return; end
    if all(ptsSelected)
        gui.myWarndlg(FigureHandle, "All cells are in the same group.");
        ptsSelected = [];
        return;
    end
else
    % FigureHandle, so an expanded selection is highlighted here too.
    [ptsSelected, letdoit] = gui.i_expandbrushed(brushed, sce, FigureHandle);
    if ~letdoit, ptsSelected = []; end
end

end


function i_Brush4TopMarkers(FigureHandle, sce, ptsSelected, uselasso)

if nargin < 4
    uselasso = true;
end

axesh = FigureHandle.findobj('type', 'Axes');
[axx, bxx] = view(axesh);

[numfig] = gui.i_inputnumg(500);
if isempty(numfig), return; end


if uselasso, fw = gui.myWaitbar(FigureHandle); end
y = double(ptsSelected);
% `sce.c = 1 + ptsSelected` was here. SCE.C is the app's active cell
% grouping and drives the colouring of the main plot; SCE is a handle
% object, so that line replaced the user's clustering or cell-type
% grouping with a two-level brushed/not-brushed indicator, permanently and
% with no way back. Nothing in this callback ever read it -- the
% regressions below use Y -- so it was pure collateral damage.
X = sce.X';

% uselasso = true;

try
    if issparse(X), X = full(X); end
%    assignin("base","X",X);
%    assignin("base","y",y);
    if uselasso
        [B] = lasso(X, y, 'DFmax', numfig*3, 'MaxIter', 1e3);
        [~, ix] = min(abs(sum(B > 0)-numfig));
        b = B(:, ix);
        idx = b > 0;
    else
        % mdl = fitglm(X, y, 'Distribution', 'binomial', 'Link', 'logit');
        % B = mdl.Coefficients.Estimate;
        % [~, idx] = mink(B, numfig);
        idx = LRDETest(X, y, numfig);
    end
catch ME
    if uselasso, gui.myWaitbar(FigureHandle, fw, true); end
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    rethrow(ME);
end

if ~any(idx)
   if uselasso, gui.myWaitbar(FigureHandle, fw); end
    gui.myWarndlg(FigureHandle, 'No marker found');
    return;
end

if uselasso
    markerlist = sce.g(idx);
    [~, jx] = sort(b(idx), 'descend');
    markerlist = markerlist(jx);
else
    markerlist = sce.g(idx);
end

fprintf('%d marker genes: ', length(markerlist));
fprintf('%s ', markerlist)
fprintf('\n')

n = length(markerlist);
y = cell(n,1);
for k=1:n
    y{k} = sce.X(sce.g == markerlist(k), :);
end

gui.sc_uitabgrpfig_expplot(y, markerlist, sce.s, FigureHandle, [axx, bxx]);
if uselasso, gui.myWaitbar(FigureHandle, fw); end


function idx = LRDETest(X, y, k)
% The likelihood-ratio ranking lives in PKG.E_LRDETEST now, so it can be
% tested without a figure. It used to be written out here with
%
%     n = size(X, 1);          % number of CELLS
%     for x = 1:size(X, 1)     % iterating over CELLS
%         model_data = table(X(:, x), ...);   % indexing a GENE
%
% and X is cells-by-genes (X = sce.X' above), so it tested genes 1..nCells
% and never looked at the rest, while line 101 below uses the returned
% indices to pick gene names. Measured on 40 cells and 120 genes with five
% planted markers at genes 100-104: not one of the five was tested, and the
% callback reported five arbitrary noise genes as the markers of the
% brushed selection. With fewer genes than cells it raised
% MATLAB:badsubscript instead.
fw = gui.myWaitbar(FigureHandle);
try
    idx = pkg.e_lrdetest(X, y, k, ...
        @(frac) gui.myWaitbar(FigureHandle, fw, false, '', '', frac));
catch ME
    gui.myWaitbar(FigureHandle, fw, true);
    rethrow(ME);
end
gui.myWaitbar(FigureHandle, fw);
end
end
