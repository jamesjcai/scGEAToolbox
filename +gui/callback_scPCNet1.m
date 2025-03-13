function callback_scPCNet1(src, events)

[FigureHandle, sce] = gui.gui_getfigsce(src);

extprogname = 'ml_scTenifoldNet';
preftagname = 'externalwrkpath';
[wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname, FigureHandle);
if isempty(wkdir), return; end


if numel(unique(sce.c_cell_type_tx)) > 1
    answer = gui.myQuestdlg(FigureHandle, 'Construct gene regulatory network (GRN) for all cells or selected cells?', ...
            '', {'All Cells', 'Select Cells...', 'Cancel'}, ...
            'All Cells');
    switch answer
        case 'Cancel'
            return;
        case 'All Cells'
        case 'Select Cells...'
            gui.callback_SelectCellsByClass(src, events);
            return;
        otherwise
            return;
    end
end

    %     answer=gui.myQuestdlg(FigureHandle, 'This analysis may take several hours. Continue?');
    %     if ~strcmpi(answer,'Yes'), return; end
    %useparallel = false;
    answer = gui.myQuestdlg(FigureHandle, 'Use parallel computing or not?', 'Parallel Computing', ...
        {'Use parallel', 'Not use parallel'}, 'Use parallel');
    switch answer
        case 'Use parallel'
            useparallel = true;
        case 'Not use parallel'
            useparallel = false;
        otherwise
            return;
    end

    try
        disp('>> [A]=sc_pcnet(sce.X);');
        X = sc_norm(sce.X);
        X = log1p(X);
        if useparallel
            fw = gui.myWaitbar(FigureHandle);
            [A] = sc_pcnetpar(X);
            gui.myWaitbar(FigureHandle, fw);
        else
            [A] = sc_pcnet(X, 3, false, false, true);
        end
    catch ME
        if useparallel
            gui.myWaitbar(FigureHandle, fw, true);
        end
        errordlg(ME.message);
        return;
    end


    try
        cd(wkdir);
        [~, tmpmat] = fileparts(tempname);
        g = sce.g;
        fprintf('Saving network (A) to %s.mat\n', tmpmat);
        save(tmpmat, 'A', 'g', '-v7.3');
    catch ME
        disp(ME.message);
    end

    % tstr=matlab.lang.makeValidName(string(datetime));
    % save(sprintf('A_%s',tstr),'A','g','-v7.3');

    if ~(ismcc || isdeployed)
        labels = {'Save network to variable named:', ...
            'Save sce.g to variable named:'};
        vars = {'A', 'g'};
        values = {A, sce.g};
        waitfor(export2wsdlg(labels, vars, values));
    end

    answer = gui.myQuestdlg(FigureHandle, 'Save network A to MAT file?');
    switch answer
        case 'Yes'
            [file, path] = uiputfile({'*.mat'; '*.*'}, 'Save as', ...
                'GeneRegulatoryNetwork');
            figure(FigureHandle);
            if isequal(file, 0) || isequal(path, 0)
                return;
            else
                filename = fullfile(path, file);
                fw = gui.myWaitbar(FigureHandle);
                g = sce.g;
                save(filename, 'A', 'g', '-v7.3');
                gui.myWaitbar(FigureHandle, fw);
            end
    end

end
