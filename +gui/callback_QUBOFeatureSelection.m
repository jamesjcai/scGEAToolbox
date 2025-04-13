function callback_QUBOFeatureSelection(src, ~)

    [FigureHandle, sce] = gui.gui_getfigsce(src);


    installedPackages = matlabshared.supportpkg.getInstalled;
    isQuantumInstalled = any(strcmp({installedPackages.Name}, 'MATLAB Support Package for Quantum Computing'));
    
    if ~isQuantumInstalled
        gui.myErrordlg(FigureHandle, 'Quantum Computing Support Package is not installed');
        return;
    end

    extprogname = 'scgeatool_QUBOFSAnalysis';
    preftagname = 'externalwrkpath';
    [wrkdir] = gui.gui_setprgmwkdir(extprogname, preftagname, FigureHandle);
    if isempty(wrkdir), return; end

    answer = gui.myQuestdlg(FigureHandle, 'Select a dependent variable y. Continue?','');
    if ~strcmp(answer,'Yes'), return; end

    [thisx, xlabelv] = gui.i_select1state(sce, false, false, false, true, FigureHandle);
    if isempty(thisx), return; end
    if ~isnumeric(thisx)
        gui.myWarndlg(FigureHandle, 'This function works with continuous varibles only.');
        return;
    end

    k = gui.i_inputnumk(20, 2, sce.NumGenes, 'Number of features (genes)', FigureHandle);
    if isempty(k), return; end

    [Xt] = gui.i_transformx(sce.X, true, 5, FigureHandle);
    if isempty(Xt), return; end
    
    gui.i_resetrngseed(src, [], false);
    fw = gui.myWaitbar(FigureHandle);
    
    b = qtm.qubofs(Xt, thisx, k);
    
    gui.myWaitbar(FigureHandle, fw);

    selected_genes = sce.g(b);
    T = table(selected_genes);
    outfile = sprintf('QUBO_Selected_%d_Genes_%s', k, xlabelv);
            [~, ~] = gui.i_exporttable(T, true, ...
                'Tqubofgenes', ...                 
                outfile, [], "Selected_Genes", FigureHandle);

        % if ~isempty(filesaved)
        %     gui.myHelpdlg(FigureHandle, sprintf('Result has been saved in %s',filesaved));
        %     %fprintf('Result has been saved in %s\n', filesaved);
        % end

end