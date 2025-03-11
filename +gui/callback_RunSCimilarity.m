function [needupdatesce] = callback_RunSCimilarity(src, ~)

    needupdatesce = false;
    [y, prepare_input_only] = gui.i_memorychecked;
    if ~y, return; end

    % [~, sys] = memory;
    % totalMemoryGB = sys.PhysicalMemory.Total / (1024^3);
    % fprintf('Total Physical Memory: %.2f GB\n', totalMemoryGB);
    % if totalMemoryGB < 32
    %     answer = gui.myQuestdlg(FigureHandle, '>=64 GB of memory is recommanded. The computer has less than 32 GB. Continue?', ...
    %         '', 'Yes, still run', 'No, prepare input only', 'Cancel', 'Yes, still run');
    %     switch answer
    %         case 'Yes, still run'
    %             prepare_input_only = false;
    %         case 'No, prepare input only'
    %             prepare_input_only = true;
    %         case 'Cancel'
    %             return;
    %         otherwise
    %             return;
    %     end
    % end
    
[FigureHandle, sce] = gui.gui_getfigsce(src);

    % https://genentech.github.io/scimilarity/notebooks/cell_annotation_tutorial.html
    % SCimilarity trained model. Download SCimilarity models. 
    % Note, this is a large tarball - downloading and uncompressing can take a several minutes.
    
    [modeldir] = gui.i_setscimilaritymodelpath;
    if isempty(modeldir), return; end

    label_ints_file = fullfile(modeldir, 'label_ints.csv');
    if exist(label_ints_file, "file")
        answer = gui.myQuestdlg(FigureHandle, 'Unconstrained or constrained annotation','', ...
            {'Unconstrained','Constrained'},'Unconstrained');
        switch answer
            case 'Unconstrained'
                target_celltypes = '';
            case 'Constrained'
                warning off
                T = readtable(label_ints_file, 'ReadVariableNames',true);
                warning on
                allcelltypes = natsort(string(T.x0));
                %scimilmodelpath
                %scimiltargetcel
                if ispref('scgeatoolbox', 'scimiltargetcel')
                    preselected_celltypes = getpref('scgeatoolbox', 'scimiltargetcel');
                else
                    preselected_celltypes = '';
                end
                [idx] = gui.i_selmultidlg(allcelltypes, preselected_celltypes, FigureHandle);
                if isempty(idx), return; end
                if idx == 0, return; end
                target_celltypes = allcelltypes(idx);
                setpref('scgeatoolbox', 'scimiltargetcel', target_celltypes);
            otherwise
                return;
        end
    else
        gui.myWarndlg(FigureHandle, "Missing label_ints.csv. Scimilarity model path is invalid.");
        return;
    end

    extprogname = 'py_scimilarity';
    preftagname = 'externalwrkpath';
    [wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname, FigureHandle);
    if isempty(wkdir), return; end
    sce.g = upper(sce.g);

    
    if prepare_input_only
        try
            fw = gui.gui_waitbar;
            run.py_scimilarity(sce, modeldir, wkdir, target_celltypes, true, prepare_input_only);
            gui.gui_waitbar(fw);
            if strcmp(gui.myQuestdlg(FigureHandle, 'Input files prepared. Open the working folder?'),'Yes')
                winopen(wkdir);
            end            
        catch ME
            gui.gui_waitbar(fw, true);
            errordlg(ME.message);
            return;
        end
        needupdatesce = false;        
    else    
        fw = gui.gui_waitbar;
        try
            [c] = run.py_scimilarity(sce, modeldir, wkdir, target_celltypes, true);
            assert(sce.NumCells==numel(c));
            if ~(isscalar(unique(sce.c_cell_type_tx)) && unique(sce.c_cell_type_tx)=="undetermined")
                sce.list_cell_attributes = [sce.list_cell_attributes, ...
                    {'old_cell_type', sce.c_cell_type_tx}];
            end
            sce.c_cell_type_tx = c;
        catch ME
            gui.gui_waitbar(fw, true);
            errordlg(ME.message);
            return;
        end    
        guidata(FigureHandle, sce);
        needupdatesce = true;
        gui.gui_waitbar(fw);
    end
end
