function [needupdate] = xcallback_RunSeuratSCTransform(src, ~)
needupdate = false;
[ok] = gui.i_confirmscript('Run Seurat sctransform to obtain normalized expression matrix?', ...
    'R_SeuratSctransform', 'r');
if ~ok, return; end

FigureHandle = src.Parent.Parent;
sce = guidata(FigureHandle);
fw = gui.gui_waitbar;
try
    [X] = run.r_SeuratSctransform(sce.X, sce.g);
catch ME
    gui.gui_waitbar(fw);
    errordlg(ME.message);
    return;
end
gui.gui_waitbar(fw);
if ~isempty(X)
    if isequal(size(X), size(sce.X))
        answer = questdlg('Update current SCE.X with transformed X or just export transformed X', '', ...
            'Update', 'Export', 'Update');
        switch answer
            case 'Update'
                needupdate = true;
                sce.X = X;
                guidata(FigureHandle, sce);
                helpdlg('SCE.X has been updated.');
            case 'Export'
                labels = {'Save transformed X to variable named:'};
                vars = {'X'};
                values = {X};
                export2wsdlg(labels, vars, values, ...
                    'Save Data to Workspace');
            otherwise
                errordlg('Invalid selection.');
        end
    end
else
    errordlg("Seurat/sctransform runtime error.");
end
end
