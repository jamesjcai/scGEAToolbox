function xcallback_TwoGeneCooccurrenceTest(~, ~)

% p=run.CooccurrenceAffinity(X)
% pkg.adjustedrandindex

    helpdlg('This function is under development.', '');
        return;

        answer = questdlg('Calculate gene expression mean, CV, and dropout rate. Save output to a table.', '');
        if ~strcmp(answer, 'Yes'), return; end
        try
            FigureHandle = src.Parent.Parent;
            sce = guidata(FigureHandle);
            fw = gui.gui_waitbar;
            [T] = sc_genestats(sce.X, sce.g)
            gui.gui_waitbar(fw);
            gui.i_exporttable(T);
        catch ME
            errordlg(ME.message);
        end
end
