function callback_BuildGeneNetwork(src, ~)
[FigureHandle, sce] = gui.gui_getfigsce(src);

[glist] = gui.i_selectngenes(sce, [], FigureHandle);
if isempty(glist), return; end

[y, i] = ismember(upper(glist), upper(sce.g));
if ~all(y), error('Runtime error.'); end
fprintf("% s\n", glist)

methods = {'PCR (PC Regression)', ...
           'Chatterjee Xi Correlation', ...
           'Pearson Correlation', ...
           'Distance Correlation', ...
           'Mutual Information', ...
           'GENIE3 (Random Forest)'};
methodkeys = {'pcrnet', 'xicor', 'pearson', 'distcorr', 'mi', 'genie3'};

[sel, ok] = gui.myListdlg(FigureHandle, methods, ...
    'Select GRN construction method', methods{1}, false);
if ~ok, return; end
methodkey = methodkeys{sel};

[Xt] = gui.i_transformx(sce.X, true, 5, FigureHandle);
if isempty(Xt), return; end
x = Xt(i, :);

fw = gui.myWaitbar(FigureHandle);
try
    A = sc_grn(x, methodkey);
catch ME
    gui.myWaitbar(FigureHandle, fw);
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    return;
end
gui.myWaitbar(FigureHandle, fw);

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
        G = net.i_makegraph(A, glist);
        if strcmp('Yes', gui.myQuestdlg(FigureHandle, 'Save network?'))
            if gui.i_isuifig(FigureHandle)
                [file, path] = uiputfile(FigureHandle, {'*.mat'; '*.*'}, ...
                    'Save as', 'network_file');
            else
                [file, path] = uiputfile({'*.mat'; '*.*'}, ...
                    'Save as', 'network_file');
            end
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
