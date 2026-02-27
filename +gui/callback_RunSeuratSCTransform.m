function [needupdate] = callback_RunSeuratSCTransform(src,~)
    needupdate=false;
    [FigureHandle, sce] = gui.gui_getfigsce(src);

    answer2 = gui.myQuestdlg(FigureHandle, ...
        'Perform SCTransform or load saved transformed X?', ...
        '', {'Perform Transform', 'Load Saved', 'Cancel'}, 'Draw Curve');
    switch answer2
        case 'Perform Transform'

        case 'Load Saved'
            if gui.i_isuifig(FigureHandle)
                [file, path] = uigetfile(FigureHandle, '*.mat', ...
                    'Select a MAT-file to Load');
            else
                [file, path] = uigetfile('*.mat', ...
                    'Select a MAT-file to Load');
            end
            if isequal(file, 0)
                disp('User canceled the file selection.');
                return;
            end
            
            fullFileName = fullfile(path, file);
            loadedData = load(fullFileName);
            if isfield(loadedData, 'X')
                X = loadedData.X;
            else
                gui.myErrordlg(FigureHandle, 'Not a valid .mat file.','');
                return;
            end
            if strcmp('Yes', gui.myQuestdlg(FigureHandle,'Transformed X has been loaded. Use it to update SCE.X?'))
               needupdate = true;
               sce.X = X;
               gui.myGuidata(FigureHandle, sce, src);
               gui.myHelpdlg(FigureHandle, 'SCE.X has been updated.');
            end
           return;            
        otherwise
           return;
    end

    extprogname = 'R_SeuratSctransform';
    preftagname = 'externalwrkpath';
    [wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname, FigureHandle);
    if isempty(wkdir), return; end
    %[ok] = gui.i_confirmscript('Run Seurat sctransform to obtain normalized expression matrix?', ...
    %    'R_SeuratSctransform','r');
    %if ~ok, return; end

    fw = gui.myWaitbar(FigureHandle);
    try
        [X, scale_X] = run.r_SeuratSctransform(sce.X, sce.g, wkdir);
    catch ME
        gui.myWaitbar(FigureHandle, fw);
        gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
        return;
    end       
    gui.myWaitbar(FigureHandle, fw);

    if ~isempty(X)
        if isequal(size(X), size(sce.X))
            answer = gui.myQuestdlg(FigureHandle, 'Update current SCE.X with transformed X or save transformed X','', ...
                {'Update','Save'}, 'Update');
            switch answer
                case 'Update'
                   needupdate = true;
                   sce.X = X;
                   gui.myGuidata(FigureHandle, sce, src);
                   gui.myHelpdlg(FigureHandle, 'SCE.X has been updated.');
                case 'Export'
                    labels = {'Save transformed X to variable named:'}; 
                    vars = {'X','scale_X'};
                    values = {X, scale_X};
                    export2wsdlg(labels,vars,values,...
                            'Save Data to Workspace');
                case 'Save'
                    if gui.i_isuifig(FigureHandle)
                        [file, path] = uiputfile(FigureHandle, '*.mat', 'Save as', ...
                            'sctransformed_X.mat');
                    else
                        [file, path] = uiputfile('*.mat', 'Save as', ...
                            'sctransformed_X.mat');
                    end
                    if isequal(file, 0) || isequal(path, 0)
                        disp('User canceled the file selection.');
                        return;
                    end                                                                        
                    fullFileName = fullfile(path, file);
                    save(fullFileName, 'X', 'scale_X');
                    disp(['Variables saved to ', fullFileName]);
                    gui.myHelpdlg(FigureHandle, sprintf('Transformed X is saved in %s.', fullFileName));
                otherwise
                    gui.myErrordlg(FigureHandle, 'Invalid selection.');
            end
        end
    else
        gui.myErrordlg(FigureHandle, "Seurat/sctransform runtime error.");
    end
    
end
