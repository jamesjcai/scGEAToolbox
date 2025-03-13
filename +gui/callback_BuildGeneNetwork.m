function callback_BuildGeneNetwork(src, ~)
    [FigureHandle, sce] = gui.gui_getfigsce(src);
    
    [glist] = gui.i_selectngenes(sce, [], FigureHandle);
    if isempty(glist), return; end
    
    [y, i] = ismember(upper(glist), upper(sce.g));
    if ~all(y), error('xxx'); end
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
            
            fw = gui.gui_waitbar_adv;
            n = size(x,1);
            A = zeros(n);
            for k=1:n
                gui.gui_waitbar_adv(fw,k/n);
                for l=1:n
                    if k~=l
                        A(k,l) = pkg.e_xicor(x(k,:),x(l,:));
                    end
                end
            end
            gui.gui_waitbar_adv(fw);
        otherwise
            return;
    end
%    [~, systemView] = memory;
%    disp(systemView.PhysicalMemory.Available)
%    bytesPerElement = 8;    % For double precision
%    maxElements = systemView.PhysicalMemory.Available / bytesPerElement;
%    maxSize = floor(sqrt(maxElements));  % If square matrix is needed
    cannotview = false;
    cannotsave = false;
    try
        sc_grnview(A, glist, [], FigureHandle);
    catch ME
        cannotview = true;
        errordlg(ME.message);
    end
    if cannotview
        try
            G = pkg.i_makegraph(A, glist);
            if strcmp('Yes', gui.myQuestdlg(FigureHandle, 'Save network?'))
                [file, path] = uiputfile({'*.mat'; '*.*'}, ...
                    'Save as', 'network_file');     
                figure(FigureHandle);
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
            errordlg(ME.message);
            cannotsave = true;
        end
    end
    if cannotview && cannotsave

    end
end
