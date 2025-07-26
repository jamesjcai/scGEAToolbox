function callback_VirtualKOGenKI(src, ~)


[FigureHandle, sce] = gui.gui_getfigsce(src);

gsorted = natsort(sce.g);
if isempty(gsorted), return; end

if gui.i_isuifig(FigureHandle)
    [indx2, tf] = gui.myListdlg(FigureHandle, gsorted, 'Select a KO gene');
else
    [indx2, tf] = listdlg('PromptString', {'Select a KO gene'}, ...
        'SelectionMode', 'single', 'ListString', ...
        gsorted, 'ListSize', [220, 300]);
end

if tf == 1
    [~, idx] = ismember(gsorted(indx2), sce.g);
else
    return;
end


answer = gui.myQuestdlg(FigureHandle, sprintf(['Ready to construct network ' ...
    'and then knock out gene #%d (%s). Continue?'], ...
    idx, sce.g(idx)));
if ~strcmpi(answer, 'Yes'), return; end

extprogname = 'py_GenKI';
preftagname = 'externalwrkpath';
[wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname, FigureHandle);
if isempty(wkdir), return; end

try
    %fw = gui.myWaitbar(FigureHandle);
    [T] = run.py_GenKI(sce.X, sce.g, idx, wkdir);
    %gui.myWaitbar(FigureHandle, fw);
catch ME
    %gui.myWaitbar(FigureHandle, fw);
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    return;
end


gui.i_exporttable(T, true, 'Tgenkiglist', 'GenKIResulTable');
% "Tcellattrib","CellAttribTable"
% "Tviolindata","ViolinPlotTable"
% "Tgenkiglist","GenKIResulTable"

disp('Downstream Analysis Options:');
disp('===============================');
disp('run.web_Enrichr(T.gene(1:200));');
disp('===============================');
end
