function callback_RunMonocle2(src,~)
    uiwait(warndlg(['Monocle 2 is obsoleted. Authors are ' ...
        'deprecating it in favor of Monocle 3.']));
    [ok]=gui.i_confirmscript('Run Monocle 2?', ...
        'R_monocle2','r');    
    if ~ok, return; end

    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);
    fw=gui.gui_waitbar;
    try
        [t_mono, s_mono] = run.r_monocle2(sce.X);
    catch ME
        gui.gui_waitbar(fw,true);
        errordlg(ME.message);
        return;
    end
    gui.gui_waitbar(fw);
    if isempty(t_mono) || isempty(s_mono)
        errordlg('MONOCLE running time error.');
        return;
    end
    if length(t_mono)~=sce.NumCells
        errordlg('MONOCLE running time error.');
        return;
    end
    [y,idx]=ismember({'monocle_pseudotime'},...
        sce.list_cell_attributes(1:2:end));
    if y   
        sce.list_cell_attributes{idx+1}=t_mono;
    else
        sce.list_cell_attributes=[sce.list_cell_attributes,...
            {'monocle_pseudotime',t_mono}];
    end
    sce.struct_cell_embeddings.('monocle')=s_mono;
    
    guidata(FigureHandle,sce);
    
    if ~(ismcc || isdeployed)    
        labels = {'Save pseudotime T to variable named:', ...
            'Save S to variable named:'};
        vars = {'t_mono', 's_mono'};
        values = {t_mono, s_mono};
        msgfig=export2wsdlg(labels, vars, values);
        uiwait(msgfig)
    else
        gui.i_exporttable(table(t_mono),true,'t_mono');
    end

    answer = questdlg('View Monocle DDRTree?', ...
        'Pseudotime View', ...
        'Yes', 'No', 'Yes');
    switch answer
        case 'Yes'
            FigureHandle=figure;
            gui.i_gscatter3(s_mono,grp2idx(t_mono));
            colorbar
            defaultToolbar = findall(FigureHandle, 'tag','FigureToolBar');  % get the figure's toolbar handle
            gui.add_3dcamera(defaultToolbar, 'Scores');
            hc = colorbar;
            hc.Label.String = 'Pseudotime';
        otherwise
            return;
    end
end
