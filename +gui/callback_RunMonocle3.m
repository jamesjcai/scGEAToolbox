function callback_RunMonocle3(src,~)

    FigureHandle=src.Parent.Parent;
    a=findall(FigureHandle,'type','axes');
    h=findall(a,'type','scatter');
    ptsSelected = logical(h.BrushData.');
    if ~any(ptsSelected)
        helpdlg('Please use the brush in the axes toolbar to select root cell(s).','');
        return;
    end
    idx=find(ptsSelected);

    [ok]=gui.i_confirmscript('Run Pseudotime Analysis (Monocle3)?', ...
        'R_monocle3','r');    
    if ~ok, return; end
    

    sce=guidata(FigureHandle);
    fw=gui.gui_waitbar;
    try
        [t_mono3] = run.r_monocle3(sce.X,idx);
    catch ME
        gui.gui_waitbar(fw,true);
        errordlg(ME.message);
        return;
    end
    gui.gui_waitbar(fw);
    if isempty(t_mono3)
        errordlg('MONOCLE3 running time error.');
        return;
    end
    if length(t_mono3)~=sce.NumCells
        errordlg('MONOCLE3 running time error.');
        return;
    end
    [y,idx]=ismember({'monocle3_pseudotime'},...
        sce.list_cell_attributes(1:2:end));
    if y
        sce.list_cell_attributes{idx+1}=t_mono3;
    else
        sce.list_cell_attributes=[sce.list_cell_attributes,...
            {'monocle3_pseudotime',t_mono3}];
    end
    %sce.struct_cell_embeddings.('monocle')=s_mono;
    
    guidata(FigureHandle,sce);
    
    if ~(ismcc || isdeployed)    
        labels = {'Save pseudotime T to variable named:', ...
            };
        vars = {'t_mono3'};
        values = {t_mono3};
        msgfig=export2wsdlg(labels, vars, values);
        uiwait(msgfig)
    else
        gui.i_exporttable(table(t_mono3),true,'t_mono3');
    end

    % answer = questdlg('View Monocle DDRTree?', ...
    %     'Pseudotime View', ...
    %     'Yes', 'No', 'Yes');
    % switch answer
    %     case 'Yes'
    %         FigureHandle=figure;
    %         gui.i_gscatter3(s_mono,grp2idx(t_mono));
    %         colorbar
    %         defaultToolbar = findall(FigureHandle, 'tag','FigureToolBar');  % get the figure's toolbar handle
    %         gui.add_3dcamera(defaultToolbar, 'Scores');
    %         hc = colorbar;
    %         hc.Label.String = 'Pseudotime';
    %     otherwise
    %         return;
    % end
end
