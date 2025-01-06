function [hFig] = i_heatscatterfig(sce, cs, posg, csname, parentfig)

%see also: gui.i_stemscatterfig

if nargin < 5, parentfig = []; end
if nargin < 4 || isempty(csname), csname = "CellScore"; end

hx=gui.myFigure;
hFig = hx.FigureHandle;

gui.i_heatscatter(sce.s, cs);
colorbar;
%cb.Label.String =  'Expression Level';


zlabel('Score value')
title(strrep(csname, '_', '\_'));

% tb = findall(hFig, 'Tag', 'FigureToolBar');
tb = uitoolbar('Parent', hFig);
% uipushtool(tb, 'Separator', 'off');
hx.addCustomButton('off', @in_saveScoreTable, "export.gif", 'Save cell score/gene expression to table');
hx.addCustomButton('on', @in_geneheatmapx, 'greenarrowicon.gif', 'Heatmap');
hx.addCustomButton('off', @in_genedotplot, 'greencircleicon.gif', 'Dot plot');
hx.addCustomButton('on', @in_viewgenenames, 'HDF_point.gif', 'Show gene names');
pkg.i_addbutton2fig(tb,'on', @in_stemplot,'icon-mat-blur-on-10.gif','Show stem plot');
%pkg.i_addbutton2fig(tb,'on',@i_viewscatter3,'icon-mat-blur-on-10.gif','Show scatter plot');
hx.show(parentfig);

    function in_stemplot(~,~)
        gui.i_stemscatterfig(sce, cs, posg, csname);
        % delete(h1);
        % h1 = gui.i_stemscatter(sce.s, cs);
    end

    % function i_viewscatter3(~, ~)
    %     figure;
    %     s = sce.s;
    %     x = s(:, 1);
    %     y = s(:, 2);
    %     if size(s, 2) >= 3
    %         z = s(:, 3);
    %         is2d = false;
    %     else
    %         z = zeros(size(x));
    %         is2d = true;
    %     end
    %     scatter3(x, y, z, 10, cs, 'filled');
    %     if is2d, view(2); end
    % end

    function in_viewgenenames(~, ~)
        [passed] = i_checkposg;
        if ~passed, return; end
        %         if isempty(posg)
        %             helpdlg('The gene set is empty. This score may not be associated with any gene set.');
        %         else
        idx = matches(sce.g, posg, 'IgnoreCase', true);
        gg = sce.g(idx);
        inputdlg(csname, ...
            '', [15, 80], ...
            {char(gg)});
        %        end
    end

    function in_saveScoreTable(~, ~)
        gui.i_exporttable(table(cs), false, ...
            char(matlab.lang.makeValidName(string(csname))));
    end
        
    function in_geneheatmapx(~, ~)
        [passed] = i_checkposg;
        if ~passed, return; end

        [thisc] = gui.i_select1class(sce);
        if ~isempty(thisc)
            gui.i_geneheatmap(sce, thisc, posg);
        end
    end
            
    function in_genedotplot(~, ~)
        [passed] = in_checkposg;
        if ~passed, return; end
        [thisc] = gui.i_select1class(sce);
        if isempty(thisc), return; end
        [c, cL] = grp2idx(thisc);
        idx = matches(posg, sce.g, 'IgnoreCase', true);
        if any(idx)
            gui.i_dotplot(sce.X, sce.g, c, cL, posg(idx));
        else
            helpdlg('No genes in this data set.','')
        end
    end
                
    function [passed] = in_checkposg
        if isempty(posg)
            passed = false;
            helpdlg('The gene set is empty. This score may not be associated with any gene set.','');
        else
            passed = true;
        end
    end

end
