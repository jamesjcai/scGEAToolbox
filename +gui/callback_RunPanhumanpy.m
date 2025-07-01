function [needupdatesce, T] = callback_RunPanhumanpy(src, ~)

    needupdatesce = false;
    %[y, prepare_input_only] = gui.i_memorychecked;
    %if ~y, return; end
    prepare_input_only = false;
   
    [FigureHandle, sce] = gui.gui_getfigsce(src);

    % https://genentech.github.io/scimilarity/notebooks/cell_annotation_tutorial.html
    % SCimilarity trained model. Download SCimilarity models. 
    % Note, this is a large tarball - downloading and uncompressing can take a several minutes.
    
    %[modeldir] = gui.i_setscimilaritymodelpath;
    %if isempty(modeldir), return; end

    %label_ints_file = fullfile(modeldir, 'label_ints.csv');
    %if exist(label_ints_file, "file")
    %    answer = gui.myQuestdlg(FigureHandle, 'Unconstrained or constrained annotation','', ...
    %        {'Unconstrained','Constrained'},'Unconstrained');
    %    switch answer
    %        case 'Unconstrained'
    %            target_celltypes = '';
    %        case 'Constrained'
    %            T = readtable(label_ints_file, 'ReadVariableNames',true, ...
    %                'VariableNamingRule', 'modify');
    %            allcelltypes = natsort(string(T.x0));
                %scimilmodelpath
                %scimiltargetcel
    %            if ispref('scgeatoolbox', 'scimiltargetcel')
    %                preselected_celltypes = getpref('scgeatoolbox', 'scimiltargetcel');
    %            else
    %                preselected_celltypes = '';
    %            end
    %            [idx] = gui.i_selmultidlg(allcelltypes, preselected_celltypes, FigureHandle);
    %            if isempty(idx), return; end
    %            if idx == 0, return; end
    %            target_celltypes = allcelltypes(idx);
    %            setpref('scgeatoolbox', 'scimiltargetcel', target_celltypes);
    %        otherwise
    %            return;
    %    end
    %else
    %    gui.myWarndlg(FigureHandle, "Missing label_ints.csv. Scimilarity model path is invalid.");
    %    return;
    %end

    extprogname = 'py_panhumanpy';
    preftagname = 'externalwrkpath';
    [wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname, FigureHandle);
    if isempty(wkdir), return; end
    sce.g = upper(sce.g);

    
    if prepare_input_only
        try
            fw = gui.myWaitbar(FigureHandle);
            run.py_panhumanpy(sce, wkdir, true, prepare_input_only);
            gui.myWaitbar(FigureHandle, fw);
            if strcmp(gui.myQuestdlg(FigureHandle, 'Input files prepared. Open the working folder?'),'Yes')
                winopen(wkdir);
            end            
        catch ME
            gui.myWaitbar(FigureHandle, fw, true);
            gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
            return;
        end
        needupdatesce = false;        
    else    
        fw = gui.myWaitbar(FigureHandle);
        try
            [c, T] = run.py_panhumanpy(sce, wkdir, true);
            assert(sce.NumCells==numel(c));
            if ~(isscalar(unique(sce.c_cell_type_tx)) && unique(sce.c_cell_type_tx)=="undetermined")
                sce.list_cell_attributes = [sce.list_cell_attributes, ...
                    {'old_cell_type', sce.c_cell_type_tx}];
            end
            sce.c_cell_type_tx = c;
        catch ME
            gui.myWaitbar(FigureHandle, fw, true);
            gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
            return;
        end    
        gui.myGuidata(FigureHandle, sce, src);
        needupdatesce = true;
        gui.myWaitbar(FigureHandle, fw);
    end
end
