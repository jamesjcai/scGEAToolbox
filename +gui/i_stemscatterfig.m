function [f0] = i_stemscatterfig(sce, cs, posg, csname)

if nargin < 4 || isempty(csname), csname = "CellScore"; end

[f0] = figure('Visible', false);


gui.i_stemscatter(sce.s, cs);

zlabel('Score value')
title(strrep(csname, '_', '\_'));


tb = uitoolbar(f0);
pkg.i_addbutton2fig(tb, 'off', @i_saveCrossTable, "export.gif", 'Save cross-table');
pkg.i_addbutton2fig(tb, 'off', {@gui.i_savemainfig, 3}, "powerpoint.gif", 'Save Figure to PowerPoint File...');
pkg.i_addbutton2fig(tb, 'on', @gui.i_pickcolormap, 'plotpicker-compass.gif', 'Pick new color map...');
pkg.i_addbutton2fig(tb, 'on', @gui.i_invertcolor, 'plotpicker-comet.gif', 'Invert colors');
pkg.i_addbutton2fig(tb, 'on', @i_geneheatmapx, 'greenarrowicon.gif', 'Heatmap');
pkg.i_addbutton2fig(tb, 'on', @i_genedotplot, 'greencircleicon.gif', 'Dot plot');
pkg.i_addbutton2fig(tb, 'on', @i_viewgenenames, 'HDF_point.gif', 'Show gene names');
%pkg.i_addbutton2fig(tb,'on',@i_viewscatter3,'icon-mat-blur-on-10.gif','Show scatter plot');

movegui(f0, 'center');
set(f0, 'Visible', true);

%     function i_viewscatter3(~,~)
%         figure;
%         s=sce.s;
%         x = s(:, 1); y = s(:, 2);
%         if size(s, 2) >= 3
%             z = s(:, 3);
%             is2d=false;
%         else
%             z = zeros(size(x));
%             is2d=true;
%         end
%         scatter3(x,y,z, 10, cs, 'filled');
%         if is2d, view(2); end
%     end

    function i_viewgenenames(~, ~)
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

        function i_saveCrossTable(~, ~)
            gui.i_exporttable(table(cs), false, ...
                char(matlab.lang.makeValidName(string(csname))));
    end
            function i_geneheatmapx(~, ~)
                [passed] = i_checkposg;
                if ~passed, return; end

                [thisc] = gui.i_select1class(sce);
                if ~isempty(thisc)
                    gui.i_geneheatmap(sce, thisc, posg);
                end
        end
                function i_genedotplot(~, ~)
                    [passed] = i_checkposg;
                    if ~passed, return; end
                    [thisc] = gui.i_select1class(sce);
                    [c, cL] = grp2idx(thisc);
                    idx = matches(posg, sce.g, 'IgnoreCase', true);
                    if any(idx)
                        gui.i_dotplot(sce.X, sce.g, c, cL, posg(idx));
                    else
                        helpdlg('No genes in this data set.')
                    end
            end
                    function [passed] = i_checkposg
                        if isempty(posg)
                            passed = false;
                            helpdlg('The gene set is empty. This score may not be associated with any gene set.');
                        else
                            passed = true;
                        end
                end
                end