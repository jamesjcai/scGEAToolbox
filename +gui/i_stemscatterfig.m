function [f0]=i_stemscatterfig(sce,cs,posg,csname)

        if nargin<4, csname="CellScore"; end

        [f0]=figure('Visible',false);
        gui.i_stemscatter(sce.s,cs);
        zlabel('Score Value')
        title(strrep(csname,'_','\_'));
        tb = uitoolbar(f0);
        pkg.i_addbutton2fig(tb,'off',@i_saveCrossTable,"export.gif",'Save cross-table');
        pkg.i_addbutton2fig(tb,'off',{@gui.i_savemainfig,3},"powerpoint.gif",'Save Figure to PowerPoint File...');
        pkg.i_addbutton2fig(tb,'on',@gui.i_pickcolormap,'plotpicker-compass.gif','Pick new color map...');
        pkg.i_addbutton2fig(tb,'on',@gui.i_invertcolor,'plotpicker-comet.gif','Invert colors');
        pkg.i_addbutton2fig(tb,'on',@i_geneheatmapx,'greenarrowicon.gif','Heatmap');
        pkg.i_addbutton2fig(tb,'on',@i_genedotplot,'greencircleicon.gif','Dot plot');
        pkg.i_addbutton2fig(tb,'on',@i_viewgenenames,'HDF_point.gif','Show gene names');

        movegui(f0,'center');
        set(f0,'Visible',true);

    function i_viewgenenames(~,~)
        idx=matches(sce.g,posg,'IgnoreCase',true);
        gg=sce.g(idx);
        inputdlg(csname, ...
            '',[10 50], ...
            {char(gg)});
    end
    function i_saveCrossTable(~,~)
        gui.i_exporttable(table(cs),false,csname);
    end
    function i_geneheatmapx(~,~)
        [thisc]=gui.i_select1class(sce);
        if ~isempty(thisc)
            gui.i_geneheatmap(sce,thisc,posg);
        end
    end
    function i_genedotplot(~,~)
        [thisc]=gui.i_select1class(sce);
        [c,cL]=grp2idx(thisc);
        idx=matches(posg,sce.g,'IgnoreCase',true);
        if any(idx)
            gui.i_dotplot(sce.X,sce.g,c,cL,posg(idx));
        else
            helpdlg('No genes in this data set.')
        end
    end
end