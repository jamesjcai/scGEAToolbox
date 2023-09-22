function [hFig] = i_crosstalkgraph(OUT, k, sce)

if nargin < 3, sce = []; end
if nargin < 2, k = 1; end

%    figure;
%    m=m-diag(diag(m));
%    h=heatmap(m);
%    h.XDisplayLabels=cL;
%    h.YDisplayLabels=cL;
% pkg.heatmap(m, cL, cL,'%0.2f', 'TickAngle', 45, 'ShowAllTicks', true, 'TextColor', 'w');
if ~isempty(sce)
    [c, cL] = grp2idx(sce.c_cell_type_tx);
    hFig = figure;
    ax1 = subplot(2, 2, 1);
    gui.i_gscatter3(sce.s(:, 1:2), c);
    title('Cell Types')
    box on
    dt = datacursormode;
    dt.UpdateFcn = {@i_myupdatefcnx};

    ax2 = subplot(2, 2, 2);
    ligand_mat = OUT.ligand_mat;
    receptor_mat = OUT.receptor_mat;
    ligandok = OUT.ligandok;
    receptorok = OUT.receptorok;
    a = ligand_mat(k, :);
    b = receptor_mat(k, :);
    m = (a' * b) .* ((a > 0)' * (b > 0));
    m = m ./ sum(m(m > 0));
    figname = sprintf('%s (ligand) -> %s (receptor)', ...
        ligandok(k), receptorok(k));
    G = pkg.i_makegraph(m, OUT.cL);
    p = plot(G);
    cc = repmat([0, 0.4470, 0.7410], G.numedges, 1);
    cc(G.Edges.Weight < 0, :) = repmat([0.8500, 0.3250, 0.0980], ...
        sum(G.Edges.Weight < 0), 1);
    p.EdgeColor = cc;
    w = 3;
    G.Edges.LWidths = abs(w*G.Edges.Weight/max(G.Edges.Weight));
    p.LineWidth = G.Edges.LWidths;
    title('Interaction Graph')
    ax3 = subplot(2, 2, 3);
    sc_scattermarker(sce.X, upper(sce.g), sce.s, OUT.ligandok(k), 1, [], false);
    box on
    ax4 = subplot(2, 2, 4);
    sc_scattermarker(sce.X, upper(sce.g), sce.s, OUT.receptorok(k), 1, [], false);
    box on

    tb = uitoolbar(hFig);
    pt5pickcolr = uipushtool(tb, 'Separator', 'off');
    [img, map] = imread(fullfile(fileparts(mfilename('fullpath')), ...
        '../', 'resources', 'fvtool_fdalinkbutton.gif')); % plotpicker-pie
    ptImage = ind2rgb(img, map);
    pt5pickcolr.CData = ptImage;
    pt5pickcolr.Tooltip = 'Link subplots';
    pt5pickcolr.ClickedCallback = @i_linksubplots;
    pkg.i_addbutton2fig(tb, 'on', {@gui.i_savemainfig, 3}, "powerpoint.gif", 'Save Figure to PowerPoint File...');
    sgtitle(figname);

    %             hFig.Position(3) = hFig.Position(3) * 1.5;
    %             hFig.Position(4) = hFig.Position(4) * 1.5;
    colormap(ax1, 'default');
else
    ligand_mat = OUT.ligand_mat;
    receptor_mat = OUT.receptor_mat;
    ligandok = OUT.ligandok;
    receptorok = OUT.receptorok;

    a = ligand_mat(k, :);
    b = receptor_mat(k, :);
    m = (a' * b) .* ((a > 0)' * (b > 0));
    m = m ./ sum(m(m > 0));

    figname = sprintf('%s (ligand) -> %s (receptor)', ...
        ligandok(k), receptorok(k));
    sc_grnview(m, OUT.cL, figname);
end


    function i_linksubplots(~, ~)
        %linkprop([ax3 ax4],{'CameraPosition','CameraUpVector'});
        evalin('base', 'h=findobj(gcf,''type'',''axes'');');
        evalin('base', 'hlink = linkprop(h(1:2),{''CameraPosition'',''CameraUpVector''});');
        evalin('base', 'rotate3d on');
end

        function [txt] = i_myupdatefcnx(~, event_obj)
            % pos = event_obj.Position;
            idx = event_obj.DataIndex;
            txt = cL(c(idx));
    end
    end
