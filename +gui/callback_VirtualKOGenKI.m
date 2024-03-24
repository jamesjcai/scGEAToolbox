function callback_VirtualKOGenKI(src, ~)

FigureHandle = src.Parent.Parent;
if ~gui.gui_showrefinfo('GenKI [PMID:37246643]'), return; end

sce = guidata(FigureHandle);

gsorted = natsort(sce.g);
if isempty(gsorted), return; end
[indx2, tf] = listdlg('PromptString', {'Select a KO gene'}, ...
    'SelectionMode', 'single', 'ListString', gsorted);
if tf == 1
    [~, idx] = ismember(gsorted(indx2), sce.g);
else
    return;
end


answer = questdlg(sprintf('Ready to construct network and then knock out gene #%d (%s). Continue?', ...
    idx, sce.g(idx)));
if ~strcmpi(answer, 'Yes'), return; end

extprogname = 'py_GenKI';
preftagname = 'externalwrkpath';
[wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname);
if isempty(wkdir), return; end

try
    %fw = gui.gui_waitbar;
    [T] = run.py_GenKI(sce.X, sce.g, idx, wkdir);
    %gui.gui_waitbar(fw);
catch ME
    %gui.gui_waitbar(fw);
    errordlg(ME.message);
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
