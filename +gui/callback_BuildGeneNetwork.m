function callback_BuildGeneNetwork(src, ~)
    [FigureHandle, sce] = gui.gui_getfigsce(src);
    
    [glist] = gui.i_selectngenes(sce, [], FigureHandle);
    if isempty(glist), return; end
    
    [y, i] = ismember(upper(glist), upper(sce.g));
    if ~all(y), error('Runtime error.'); end
    fprintf("%s\n", glist)
    

    switch gui.myQuestdlg(FigureHandle, 'Select algorithm:','',...
            {'PC Regression','Chaterjee Correlation'},'PC Regression')
        case 'PC Regression'
            [Xt] = gui.i_transformx(sce.X, true, 5, FigureHandle);
            if isempty(Xt), return; end
            
            x = Xt(i, :);
            
            fw = gui.myWaitbar(FigureHandle);
            A = sc_pcnet(x);
            gui.myWaitbar(FigureHandle, fw);
        case 'Chaterjee Correlation'
            [Xt] = gui.i_transformx(sce.X, true, 3, FigureHandle);
            if isempty(Xt), return; end
            
            x = Xt(i, :);
            
            fw = gui.myWaitbar(FigureHandle);
            n = size(x,1);
            A = zeros(n);
            for k=1:n
                gui.myWaitbar(FigureHandle, fw, false, '', '', k/n);
                for l=1:n
                    if k~=l
                        A(k,l) = pkg.e_xicor(x(k,:),x(l,:));
                    end
                end
            end
            gui.myWaitbar(FigureHandle, fw);
        otherwise
            return;
    end

    cannotview = false;
    cannotsave = false;
    try
        sc_grnview(A, glist, '', FigureHandle);
    catch ME
        cannotview = true;
        gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    end
    if cannotview
        try
            G = pkg.i_makegraph(A, glist);
            if strcmp('Yes', gui.myQuestdlg(FigureHandle, 'Save network?'))
                [file, path] = uiputfile({'*.mat'; '*.*'}, ...
                    'Save as', 'network_file');     
                if isvalid(FigureHandle) && isa(FigureHandle, 'matlab.ui.Figure'), figure(FigureHandle); end
                if isequal(file, 0) || isequal(path, 0)
                    return;
                else
                    filename = fullfile(path, file);
                    fw = gui.myWaitbar(FigureHandle);
                    save(filename, 'G');
                    gui.myWaitbar(FigureHandle, fw);
                    gui.myHelpdlg(FigureHandle, 'File saved.');                    
                end
            end
        catch ME
            gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
            cannotsave = true;
        end
    end
    if cannotview && cannotsave

    end
end
