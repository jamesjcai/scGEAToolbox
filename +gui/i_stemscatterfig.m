function [hFig] = i_stemscatterfig(sce, cs, posg, csname, parentfig)

%see also: gui.i_heatscatterfig

if nargin < 5, parentfig = []; end
if nargin < 4 || isempty(csname), csname = "CellScore"; end

hx = gui.myFigure(parentfig);
hFig = hx.FigHandle;
ax = hx.AxHandle;
if isempty(ax), ax = gca; end

gui.i_stemscatter(sce.s, cs, ax);

zlabel(ax, 'Score value');
title(ax, strrep(csname, '_', '\_'));

hx.addCustomButton('off', @i_saveCrossTable, "floppy-disk-arrow-in.jpg", 'Save cross-table');
hx.addCustomButton('on', @in_geneheatmapx, 'greenarrowicon.gif', 'Heatmap');
hx.addCustomButton('on', @i_genedotplot, 'greencircleicon.gif', 'Dot plot');
hx.addCustomButton('on', @i_viewgenenames, 'HDF_point.gif', 'Show gene names');
hx.addCustomButton('on', @in_heatscatterplot,'icon-mat-blur-on-10.gif','Show heat scatter plot');
hx.show(parentfig);

    function in_heatscatterplot(~, ~)
        gui.i_heatscatterfig(sce, cs, posg, csname, hFig);
        % delete(h1);
        % h1 = gui.i_stemscatter(sce.s, cs);
    end


    function i_viewgenenames(~, ~)
        [passed] = i_checkposg;
        if ~passed, return; end

        idx = matches(sce.g, posg, 'IgnoreCase', true);
        gg = sce.g(idx);
        if gui.i_isuifig(parentfig)
            answer = gui.myInputdlg({csname}, '', {char(gg)}, parentfig);
        else
            answer = inputdlg(csname, '', [15, 80], {char(gg)});
        end
        %        end
    end

    function i_saveCrossTable(~, ~)
        gui.i_exporttable(table(cs), false, ...
            char(matlab.lang.makeValidName(string(csname))));
    end

    function in_geneheatmapx(~, ~)
        [passed] = i_checkposg;
        if ~passed, return; end

        [thisc] = gui.i_select1class(sce,[],[],[],parentfig);
        if ~isempty(thisc)
            gui.i_geneheatmap(sce, thisc, posg);
        end
    end

    function i_genedotplot(~, ~)
        [passed] = i_checkposg;
        if ~passed, return; end
        [thisc] = gui.i_select1class(sce,[],[],[],parentfig);
        if isempty(thisc), return; end
        [c, cL] = grp2idx(thisc);
        idx = matches(posg, sce.g, 'IgnoreCase', true);
        if any(idx)
            gui.i_dotplot(sce.X, sce.g, c, cL, posg(idx));
        else
            gui.myHelpdlg(parentfig, 'No genes in this data set.')
        end
    end

    function [passed] = i_checkposg
        if isempty(posg)
            passed = false;
            gui.myHelpdlg(parentfig, ['The gene set is empty. This score ' ...
                'may not be associated with any gene set.']);
        else
            passed = true;
        end
    end

end
