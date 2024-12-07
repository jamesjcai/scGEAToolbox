function [needupdatesce] = callback_RunScimilarity(src, ~)


    needupdatesce = false;
    [~, sys] = memory;
    totalMemoryGB = sys.PhysicalMemory.Total / (1024^3);
    fprintf('Total Physical Memory: %.2f GB\n', totalMemoryGB);
    if totalMemoryGB < 32
        answer = questdlg('>=64 GB of memory is recommanded. The computer has less than 32 GB. Continue?','');
        if ~strcmp(answer, 'Yes'), return; end
    end   
    
    FigureHandle = src.Parent.Parent;
    sce = guidata(FigureHandle);

    [modeldir] = gui.i_setscimilaritymodelpath;
    if isempty(modeldir), return; end

    label_ints_file = fullfile(modeldir, 'label_ints.csv');
    if exist(label_ints_file, "file")        
        answer = questdlg('Unconstrained or constrained annotation','', ...
            'Unconstrained','Constrained','Unconstrained');
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
        warndlg("Missing label_ints.csv. Scimilarity model path is invalid.");
        return;
    end


    extprogname = 'py_scimilarity';
    preftagname = 'externalwrkpath';
    [wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname);
    if isempty(wkdir), return; end

    fw = gui.gui_waitbar;

    
    try
        sce.g = upper(sce.g);
        [c] = run.py_scimilarity(sce, modeldir, wkdir, target_celltypes);
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
    gui.gui_waitbar(fw);
    guidata(FigureHandle, sce);
    needupdatesce = true;
end


