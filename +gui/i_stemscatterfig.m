function i_stemscatterfig(sce,cs,csname)

if nargin<3, csname="CellScore"; end
        f0=figure('Visible',false);
        gui.i_stemscatter(sce.s,cs);
        zlabel('Score Value')
        title(strrep(csname,'_','\_'));
        tb = uitoolbar(f0);
        pkg.i_addbutton2fig(tb,'off',@i_saveCrossTable,"export.gif",'Save cross-table');
        pkg.i_addbutton2fig(tb,'off',{@gui.i_savemainfig,3},"powerpoint.gif",'Save Figure to PowerPoint File...');
        pkg.i_addbutton2fig(tb,'on',@gui.i_pickcolormap,'plotpicker-compass.gif','Pick new color map...');
        pkg.i_addbutton2fig(tb,'on',@gui.i_invertcolor,'plotpicker-comet.gif','Invert colors');
        pkg.i_addbutton2fig(tb,'on',@i_geneheatmapx,'plotpicker-cometx.gif','Heatmap');
        movegui(f0,'center');
        set(f0,'Visible',true);

    function i_saveCrossTable(~,~)
        gui.i_exporttable(table(cs),false,csname);
    end
    function i_geneheatmapx(~,~)
        gui.i_geneheatmap(sce,sce.c_cell_type_tx,posg);
    end
end