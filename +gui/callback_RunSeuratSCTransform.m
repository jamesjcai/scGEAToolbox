function [needupdate] = callback_RunSeuratSCTransform(src,~)
%CALLBACK_RUNSEURATSCTRANSFORM Run SCTransform v2 and hand back the result.
%
%   SCE.X holds raw counts, and every analysis normalizes on its own terms
%   when it needs to (GUI.I_TRANSFORMX). So the two matrices this produces
%   are saved or sent to the workspace rather than written back into the
%   object: X is the corrected counts on the log1p scale, matching Seurat's
%   SCT "data" slot, and scale_X is the Pearson residual matrix, which is
%   what downstream PCA and clustering use.
%
%   NEEDUPDATE is always false - nothing here changes SCE - and is returned
%   only so the signature matches the other GUI.CALLBACK_* functions.
%
%   See also GUI.I_TRANSFORMX, SC_SCTRANSFORMV2.

needupdate = false;
[FigureHandle, sce] = gui.gui_getfigsce(src);

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
        [scale_X, ~, Xcorrected] = sc_sctransformv2(sce.X);
        X = log1p(Xcorrected);
    end
catch ME
    gui.myWaitbar(FigureHandle, fw, true);
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    return;
end
gui.myWaitbar(FigureHandle, fw);

if isempty(X)
    gui.myErrordlg(FigureHandle, "Seurat/sctransform runtime error.");
    return;
end

answer = gui.myQuestdlg(FigureHandle, ...
    ['Transformed matrices are ready. SCE.X keeps the raw counts, ' ...
    'so save them to a file or send them to the workspace.'], '', ...
    {'Save to File', 'Send to Workspace', 'Cancel'}, 'Save to File');
switch answer
    case 'Save to File'
        [file, path] = uiputfile('*.mat', 'Save as', 'sctransformed_X.mat');
        if isequal(file, 0) || isequal(path, 0)
            disp('User canceled the file selection.');
            return;
        end
        fullFileName = fullfile(path, file);
        save(fullFileName, 'X', 'scale_X');
        disp(['Variables saved to ', fullFileName]);
        gui.myHelpdlg(FigureHandle, ...
            sprintf('Transformed X is saved in %s.', fullFileName));
    case 'Send to Workspace'
        labels = {'Corrected counts, log1p scale:', 'Pearson residuals:'};
        vars = {'X', 'scale_X'};
        values = {X, scale_X};
        export2wsdlg(labels, vars, values, 'Save Data to Workspace');
    otherwise
        % Cancel, or the dialog was dismissed.
        return;
end

end
