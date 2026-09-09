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

% Native MATLAB by default. SC_SCTRANSFORMV2 reproduces Seurat's v2 residuals
% to a per-gene correlation of 0.999996 and runs about ten times faster, so
% R is offered only when it is already configured.
useR = false;
hasR = ispref('scgeatoolbox', 'rexecutablepath') && ...
    ~isempty(getpref('scgeatoolbox', 'rexecutablepath', []));
if hasR
    backend = gui.myQuestdlg(FigureHandle, 'Choose SCTransform backend:', ...
        '', {'MATLAB (native)', 'R (Seurat)'}, 'MATLAB (native)');
    if isempty(backend), return; end
    useR = strcmp(backend, 'R (Seurat)');
end

if useR
    extprogname = 'R_SeuratSctransform';
    preftagname = 'externalwrkpath';
    [wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname, FigureHandle);
    if isempty(wkdir), return; end
end

fw = gui.myWaitbar(FigureHandle);
try
    if useR
        [X, scale_X] = run.r_SeuratSctransform(sce.X, sce.g, wkdir);
    else
        % scale_X is the Pearson residual matrix, which is what downstream
        % PCA and clustering use. X is the corrected counts on the log1p
        % scale, matching Seurat's SCT "data" slot, so that whatever
        % replaces SCE.X stays non-negative and count-like.
        [scale_X, ~, Xcorrected] = sc_sctransformv2(sce.X);
        X = log1p(Xcorrected);
    end
catch ME
    gui.myWaitbar(FigureHandle, fw, true);
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
                % if gui.i_isuifig(FigureHandle)
                %     [file, path] = uiputfile(FigureHandle, '*.mat', 'Save as');
                % else
                    [file, path] = uiputfile('*.mat', 'Save as', ...
                        'sctransformed_X.mat');
                %end
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
