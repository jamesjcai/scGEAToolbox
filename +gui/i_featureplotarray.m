function [hFig] = i_featureplotarray(X, g, c, cL, tgene, uselog, titletxt, parentfig)

if nargin < 8, parentfig = []; end
if nargin < 7, titletxt = []; end
if nargin < 6, uselog = false; end

cL = strrep(cL(:),'_','\_');
[yes] = ismember(tgene, g);
if ~any(yes), warning('No genes found.'); return; end
z = length(tgene) - sum(yes);
if z > 0, fprintf('%d gene(s) not in the list are excluded.\n', z); end
tgene = tgene(yes);
if issparse(X), X = full(X); end
if uselog, X = log1p(X); end

xgroupdata = categorical(cL(c));
hFig = figure('Visible', 'off', 'DockControls', 'off');
for kx = 1:length(tgene)
    nexttile;
    ydata = X(g == tgene(kx), :);
    violinplot(xgroupdata, ydata);
    box on
    title(tgene(kx));
end
if ~isempty(titletxt)
    titletxt = strrep(titletxt,'_','\_');
    sgtitle(titletxt);
end

gui.i_movegui2parent(hFig, parentfig);
if nargout < 1
    set(hFig, 'visible', 'on');
end
