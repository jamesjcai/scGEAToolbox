function i_markergenespanel(X, genelist, s, markerlist, numfig, N, ax, bx, sptltxt, parentfig)
% I_MARKERGENESPANEL Expression of each marker on the embedding, N per page.
%   One window, one tab per page of N (9 or 16) genes. AX and BX are the
%   view azimuth and elevation to apply to 3-D plots.

if nargin < 10, parentfig = []; end
if nargin < 9, sptltxt = ''; end
if nargin < 8, bx = []; end
if nargin < 7, ax = []; end
if nargin < 6, N = 9; end
n = length(markerlist);
if nargin < 5 || isempty(numfig)
    numfig = ceil(n/N);
end
if ~ismember(N, [9, 16])
    error('gui:i_markergenespanel:badPageSize', ...
        'Genes per page must be 9 or 16, not %d.', N);
end
if numfig < 1
    error('gui:i_markergenespanel:noGenes', 'There are no genes to show.');
end
% Ten pages at most. More used to raise a bare error('error'); now the
% first ten pages are drawn and the rest are named in a warning.
maxpages = 10;
if numfig > maxpages
    warning('gui:i_markergenespanel:tooManyGenes', ...
        'Showing the first %d of %d genes (%d pages of %d).', ...
        maxpages*N, n, maxpages, N);
    numfig = maxpages;
end

% Pages were separate windows, cascaded over the main app.
hx = gui.myFigure(parentfig);
delete(hx.AxHandle);
tabgp = uitabgroup(hx.FigHandle);
for kkk = 1:numfig
    tab = uitab(tabgp, 'Title', sprintf('Page %d', kkk));
    tl = tiledlayout(tab, sqrt(N), sqrt(N));
    for kk = 1:min([N, n - (kkk - 1) * N])
        in_drawgene(nexttile(tl), markerlist(kk+N*(kkk - 1)));
    end
    if ~isempty(sptltxt)
        title(tl, sptltxt);
    end
end
hx.show(parentfig);

    function in_drawgene(h1, targetg)
        % SC_SCATTERMARKER's method 2, drawn into H1 rather than gca.
        c = full(X(strcmp(genelist, targetg), :));
        if size(s, 2) > 2
            scatter3(h1, s(:, 1), s(:, 2), s(:, 3), 5, c, 'filled');
        else
            scatter(h1, s(:, 1), s(:, 2), 5, c, 'filled');
        end
        set(h1, 'XTickLabel', [], 'YTickLabel', [], 'ZTickLabel', []);
        grid(h1, "on");
        title(h1, targetg);
        subtitle(h1, gui.i_getsubtitle(c));
        if ~isempty(ax) && ~isempty(bx)
            view(h1, ax, bx);
        end
    end
end
