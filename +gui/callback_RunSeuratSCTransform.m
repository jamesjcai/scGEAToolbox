function [needupdate] = callback_RunSeuratSCTransform(src,~)
    needupdate=false;
    [FigureHandle] = gui.gui_getfigsce(src);

    extprogname = 'R_SeuratSctransform';
    preftagname = 'externalwrkpath';
    [wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname, FigureHandle);
    if isempty(wkdir), return; end
    [ok] = gui.i_confirmscript('Run Seurat sctransform to obtain normalized expression matrix?', ...
        'R_SeuratSctransform','r');
    if ~ok, return; end

    sce = guidata(FigureHandle);
    fw = gui.myWaitbar(FigureHandle);
    try
        [X] = run.r_SeuratSctransform(sce.X, sce.g, wkdir);
    catch ME
        gui.myWaitbar(FigureHandle, fw);
        gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
        return;
    end       
    gui.myWaitbar(FigureHandle, fw);

    if ~isempty(X)
        if isequal(size(X), size(sce.X))
            answer=questdlg('Update current SCE.X with transformed X or just export transformed X','', ...
                'Update','Export','Update');
            switch answer
                case 'Update'
                   needupdate = true;
                   sce.X = X;
                   guidata(FigureHandle, sce);
                   gui.myHelpdlg(FigureHandle, 'SCE.X has been updated.');
                case 'Export'
                    labels = {'Save transformed X to variable named:'}; 
                    vars = {'X'};
                    values = {X};
                    export2wsdlg(labels,vars,values,...
                            'Save Data to Workspace');
                otherwise                    
                    gui.myErrordlg(FigureHandle, 'Invalid selection.');
            end
        end
    else
        gui.myErrordlg(FigureHandle, "Seurat/sctransform runtime error.");
    end
end
