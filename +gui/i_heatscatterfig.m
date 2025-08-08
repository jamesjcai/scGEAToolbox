function [hFig] = i_heatscatterfig(sce, cs, posg, csname, parentfig)

%see also: gui.i_stemscatterfig

if nargin < 5, parentfig = []; end
if nargin < 4 || isempty(csname), csname = "CellScore"; end

hx=gui.myFigure;
hFig = hx.FigHandle;

gui.i_heatscatter(sce.s, cs, hx.AxHandle);
colorbar(hx.AxHandle);
%cb.Label.String =  'Expression Level';


zlabel(hx.AxHandle, 'Score value')
title(hx.AxHandle, strrep(csname, '_', '\_'));
hx.addCustomButton('off', @in_callback_saveScoreTable, "floppy-disk-arrow-in.jpg", 'Save cell score/gene expression to table');
hx.addCustomButton('on', @in_callback_geneheatmapx, 'greenarrowicon.gif', 'Heatmap');
hx.addCustomButton('off', @in_callback_genedotplot, 'greencircleicon.gif', 'Dot plot');
hx.addCustomButton('on', @in_callback_viewgenenames, 'HDF_point.gif', 'Show gene names');
hx.addCustomButton('on', @in_callback_stemplot,'icon-mat-blur-on-10.gif','Show stem plot');
hx.show(parentfig);

    function in_callback_stemplot(~,~)
        gui.i_stemscatterfig(sce, cs, posg, csname);
        % delete(h1);
        % h1 = gui.i_stemscatter(sce.s, cs);
    end

    function in_callback_viewgenenames(~, ~)
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

    function in_callback_saveScoreTable(~, ~)
        gui.i_exporttable(table(cs), false, ...
            char(matlab.lang.makeValidName(string(csname))));
    end
        
    function in_callback_geneheatmapx(~, ~)
        [passed] = i_checkposg;
        if ~passed, return; end

        [thisc] = gui.i_select1class(sce,[],[],[],parentfig);
        if ~isempty(thisc)
            gui.i_geneheatmap(sce, thisc, posg);
        end
    end
            
    function in_callback_genedotplot(~, ~)
        [passed] = in_callback_checkposg;
        if ~passed, return; end
        [thisc] = gui.i_select1class(sce,[],[],[],parentfig);
        if isempty(thisc), return; end
        [c, cL] = findgroups(string(thisc));
        idx = matches(posg, sce.g, 'IgnoreCase', true);
        if any(idx)
            gui.i_dotplot(sce.X, sce.g, c, cL, posg(idx));
        else
            gui.myHelpdlg(parentfig, 'No genes in this data set.');
        end
    end
                
    function [passed] = in_callback_checkposg
        if isempty(posg)
            passed = false;
            gui.myHelpdlg(parentfig, ...
                ['The gene set is empty. This score may not' ...
                ' be associated with any gene set.']);
        else
            passed = true;
        end
    end

end
