function callback_RevelioCellCycle(src,~)
    FigureHandle=src.Parent.Parent;    
    sce=guidata(FigureHandle);
        
        fw=gui.gui_waitbar;
        try
            [dc,T]=run.Revelio(sce);
            if isempty(dc)
                gui.gui_waitbar(fw);
                errordlg("Revelio Running Error");
                return;
            end
        catch ME
            gui.gui_waitbar(fw);
            errordlg(ME.message);
            rethrow(ME);
        end 
            gui.gui_waitbar(fw);

        figure;
        gscatter(dc(:,1),dc(:,2),T.Var1);        
        %defaultToolbar = findall(FigureHandle, 'tag','FigureToolBar');  % get the figure's toolbar handle
        %gui.add_3dcamera(defaultToolbar, 'RevelioDc');
        
    if ~(ismcc || isdeployed)
        labels = {'Save score values to variable named:','Save score table to variable named:'};
        vars = {'RevelioDc','RevelioTable'};
        values = {dc,T};
        export2wsdlg(labels,vars,values);
    else
        gui.i_exporttable(T,false,'RevelioTable');
    end
end
