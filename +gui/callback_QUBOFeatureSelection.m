function callback_QUBOFeatureSelection(src, ~)

[FigureHandle, sce] = gui.gui_getfigsce(src);

if ~(ismcc || isdeployed)
    %#exclude matlabshared.supportpkg.getInstalled
    installedPackages = matlabshared.supportpkg.getInstalled;
    isQuantumInstalled = any(strcmp({installedPackages.Name}, 'MATLAB Support Package for Quantum Computing'));
else
    isQuantumInstalled = false; % Default to false if not installed
end

if ~isQuantumInstalled
    gui.myErrordlg(FigureHandle, 'Quantum Computing Support Package is not installed');
    return;
end

% No working folder: nothing here writes one. It used to be asked for all
% the same, and answering Yes to its "Overwrite?" deleted the files in it.

answer = gui.myQuestdlg(FigureHandle, 'Select a dependent variable y. Continue?','');
if ~strcmp(answer,'Yes'), return; end

[thisx, xlabelv] = gui.i_select1state(sce, false, false, false, true, FigureHandle);
if isempty(thisx), return; end
if ~isnumeric(thisx)
    gui.myWarndlg(FigureHandle, 'This function works with continuous variables only.');
    return;
end

k = gui.i_inputnumk(20, 2, sce.NumGenes, ...
'Number of features (genes)', FigureHandle);
if isempty(k), return; end

[Xt] = gui.i_transformx(sce.X, true, "pearson_residuals", FigureHandle);
if isempty(Xt), return; end

if ~gui.i_resetrngseed(src, [], false), return; end   % cancelled
fw = gui.myWaitbar(FigureHandle);
try
    b = qtm.qubofs(Xt, thisx, k);
catch ME
    % The bar used to stay open over the app after a failure.
    gui.myWaitbar(FigureHandle, fw, true);
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    return;
end
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
