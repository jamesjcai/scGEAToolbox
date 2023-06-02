function callback_MELDPerturbationScore(src,~)

uiwait(helpdlg('The goal of MELD is to identify populations of cells that are most affected by an experimental perturbation. Rather than clustering the data first and calculating differential abundance of samples within clusters, MELD provides a density estimate for each scRNA-seq sample for every cell in each dataset. Comparing the ratio between the density of each sample provides a quantitative estimate the effect of a perturbation at the single-cell level. We can then identify the cells most or least affected by the perturbation.', ...
    'MELD [PMID:33558698]'));
[ok]=gui.i_confirmscript('Run MELD Perturbation Score (MELD)?', ...
    'py_MELD','python');
if ~ok, return; end

    FigureHandle=src.Parent.Parent;    
    sce=guidata(FigureHandle);
    if numel(unique(sce.c_batch_id))<2
        warndlg('No batch effect (SCE.C_BATCH_ID is empty)');
        return;
    end
    
    if ~gui.i_setpyenv
        return; 
    end
        
        
        fw=gui.gui_waitbar;
        try
            [score,T]=run.MELD(sce.X,sce.c_batch_id);
            if isempty(score) || size(score,1)~=size(sce.X,2)
                gui.gui_waitbar(fw);
                errordlg("MELD Running Error");
                return;
            end
        catch ME
            gui.gui_waitbar(fw);
            errordlg(ME.message);
            rethrow(ME);
        end 
            gui.gui_waitbar(fw);

        FigureHandle=figure;
        gui.i_gscatter3(sce.s,score(:,2));
        colorbar
        defaultToolbar = findall(FigureHandle, 'tag','FigureToolBar');  % get the figure's toolbar handle
        gui.add_3dcamera(defaultToolbar, 'MELD_Scores');
        
        if ~(ismcc || isdeployed)
            labels = {'Save score values to variable named:','Save score table to variable named:'};
            vars = {'MELDScores','MELDTable'};
            values = {score,T};
            export2wsdlg(labels,vars,values);
        else
            gui.i_exporttable(T,false,'MELDTable');
        end

end
