function i_markergeneshtml(sce, markerlist, numfig, axbx, nametag, pselected)
if nargin < 6, pselected = []; end
if nargin < 5, nametag = ''; end
if nargin < 4, axbx = []; end
if nargin < 3, numfig = 20; end
% markerlist=markerlist{1};

numfig = min([numfig, length(markerlist)]);
dirtxt = tempdir;
% dirtxt='./';

a = autumn;
a(1, :) = [.8, .8, .8];

if isempty(nametag)
    htmlfilenamex = tempname;
else
    nametag = matlab.lang.makeValidName(nametag);
    htmlfilenamex = [dirtxt, nametag];
end

fname = [htmlfilenamex, '.html'];
fid = fopen(fname, 'w');
htmlstr = "";
for k = 1:numfig
    targeetg = markerlist(k);
    h = figure('Visible', 'off');
    %         sc_scattermarker(sce.X,sce.g,sce.s,...
    %             markerlist(k),3,5,false);
    c = log2(1+sce.X(sce.g == targeetg, :));

    x = sce.s(:, 1);
    y = sce.s(:, 2);
    z = sce.s(:, 3);

    %subplot(2,2,1)
    scatter3(x, y, z, 5, c, 'filled');
    if ~isempty(pselected)
        hold on
        scatter3(sce.s(pselected, 1), sce.s(pselected, 2), ...
            sce.s(pselected, 3), 10, [.5, .5, .5]);
    end
    colormap(h, a);
    colorbar;
    title(targeetg)
    if ~isempty(axbx)
        view(axbx(1), axbx(2));
    end

    %subplot(2,2,2);
    h2 = figure('Visible', 'off');
    %             stem3(x,y,c,'marker','none','color','m');
    %             hold on
    %             scatter3(x,y,zeros(size(y)),5,c,'filled');
    scatter3(x, y, z, 5, c, 'filled');
    colormap(h2, a);
    colorbar;
    title(targeetg)
    if ~isempty(axbx)
        view(axbx(1), axbx(2));
    end


    imgfname1 = sprintf('heatmap_%s.png', targeetg);
    saveas(h, sprintf('%s%s', dirtxt, imgfname1));
    close(h);

    imgfname2 = sprintf('heatmap2%s.png', targeetg);
    saveas(h2, sprintf('%s%s', dirtxt, imgfname2));
    close(h2);

    h = figure('Visible', 'off');
    pkg.i_violinplot_groupordered(c, sce.c, ["1", "2"]);
    ylabel('log2(UMI+1)');
    title(targeetg)
    xtickangle(-45);
    imgfname3 = sprintf('violin_%s.png', targeetg);
    saveas(h, sprintf('%s%s', dirtxt, imgfname3));
    close(h);

    aax = sprintf('<center><table>\n<tr>\n<th><img src="%s" height=250></th>\n<th><img src="%s" height=250></th>\n<th><img src="%s" height=250></th>\n</tr>\n</table></center>\n', ...
        imgfname2, imgfname1, imgfname3);
    htmlstr = sprintf('%s\n%s', htmlstr, aax);
end
fprintf(fid, '<pre>\n');
fprintf(fid, '%s, ', markerlist);
fprintf(fid, '</pre>\n');
fprintf(fid, '%s', htmlstr);
fclose(fid);
pause(1);
% winopen(dirtxt);
web(fname, '-browser');
end
