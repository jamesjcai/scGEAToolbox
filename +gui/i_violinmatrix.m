function [hFig] = i_violinmatrix(X, g, c, cL, tgene, uselog, titletxt, parentfig)

% see also: i_violinarray
if nargin < 8, parentfig = []; end
if nargin < 7, titletxt = []; end
if nargin < 6, uselog = false; end

[yes] = ismember(tgene, g);
if ~any(yes), return; end
z = length(tgene) - sum(yes);
if z > 0
    fprintf('%d gene(s) not in the list are excluded.\n', z);
end
tgene = tgene(yes);

M = length(tgene);
% N=length(cL);

if issparse(X), X = full(X); end
if uselog, X = log1p(X); end
cL = strrep(cL(:),'_','\_');
% xgroupdata = categorical(cL(c));
xgroupdata = cL(c);
hFig = figure('Visible', 'off', 'DockControls', 'off');
tb = uitoolbar('Parent', hFig);
pkg.i_addbutton2fig(tb, 'on', {@gui.i_resizewin, hFig}, 'HDF_pointx.gif', 'Resize Plot Window');


for k = 1:M
    %gui.gui_waitbar_adv(fw, k / M);
    ydata = X(g == tgene(k), :);
    nexttile;   % subplot(M, 1, ct)
    pkg.violinplot(ydata.', xgroupdata, 'showdata', false, 'GroupOrder', cL);
    box on;
    title(tgene(k));    
end
if ~isempty(titletxt)
    titletxt = strrep(titletxt,'_','\_');
    sgtitle(titletxt);
end
gui.i_movegui2parent(hFig, parentfig);
if nargout < 1
    set(hFig, 'visible', 'on');
end