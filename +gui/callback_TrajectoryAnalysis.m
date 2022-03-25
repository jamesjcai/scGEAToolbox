function callback_TrajectoryAnalysis(src,~)
        answer = questdlg('Run pseudotime analysis (Monocle)?','', ...
            'Yes','Review R Script','Cancel','Yes');
        switch answer
            case 'Cancel'
                return;
            case 'Yes'

            case 'Review R Script'
                folder=fileparts(mfilename('fullpath'));
                scriptfile=fullfile(folder,'..','+run','external', ...
                    'R_monocle','script.R')
                t=fileread(scriptfile);
                %LF=char(10); 
                CR=char(13);  %  carriage return character equivalent to char(13) or sprintf('\r').
                t=strrep(t,[CR newline],newline);
                inputdlg('Review Script:','R Code',[10 90],{t});
                return;
            otherwise
                return;
        end

    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);
    fw=gui.gui_waitbar;
    [t_mono, s_mono] = run.monocle(sce.X);
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
            gui.add_3dcamera(defaultToolbar, 'MELD_Scores');
            hc = colorbar;
            hc.Label.String = 'Pseudotime';
        otherwise
            return;
    end
end