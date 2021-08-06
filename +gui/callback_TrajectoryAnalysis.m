function callback_TrajectoryAnalysis(src,~)
        answer = questdlg('Run pseudotime analysis (Monocle)?');
        if ~strcmp(answer, 'Yes')
            return
        end

    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);
    fw=gui.gui_waitbar;
    [t_mono, s_mono] = run.monocle(sce.X);
    cs=pkg.e_cellscores(sce.X,sce.g,"T_Cell_Exhaustion");
    gui.gui_waitbar(fw);

    
    labels = {'Save pseudotime T to variable named:', ...
        'Save S to variable named:'};
    vars = {'t_mono', 's_mono'};
    values = {t_mono, s_mono};
    msgfig=export2wsdlg(labels, vars, values);
    uiwait(msgfig)
    answer = questdlg('View Monocle DDRTree?', ...
        'Pseudotime View', ...
        'Yes', 'No', 'Yes');
    switch answer
        case 'Yes'
            FigureHandle=figure;
            gui.i_gscatter3(s_mono,grp2idx(t_mono));
            colorbar
            defaultToolbar = findall(FigureHandle, 'tag','FigureToolBar');  % get the figure's toolbar handle
            gui.add_3dcamera(defaultToolbar, 'MELD_Scores');
            hc = colorbar;
            hc.Label.String = 'Pseudotime';
        otherwise
            return;
    end
end