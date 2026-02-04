function callback_RunCellBender(src, ~)

[FigureHandle, ~] = gui.gui_getfigsce(src);
if ~gui.gui_showrefinfo('CellBender [PMID:37550580]', FigureHandle), return; end

extprogname = 'py_cellbender';
preftagname = 'externalwrkpath';
[wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname, FigureHandle);
if isempty(wkdir), return; end

answer=gui.myQuestdlg(FigureHandle, 'Select raw/unfiltered 10x Genomics H5 file (raw_feature_bc_matrix.h5)?', '');
if ~strcmp(answer,'Yes'), return; end

[filenm, pathname] = uigetfile( ...
    {'*.h5;*.hdf5', 'HDF5 Files (*.h5)'; ...
    '*.*', 'All Files (*.*)'}, ...
    'Pick 10x Genomics H5 file(s)','MultiSelect','off');
if isequal(filenm, 0), return; end
input_h5 = fullfile(pathname, filenm);
if ~exist(input_h5,"file")
    return;
end

% [ok] = gui.i_confirmscript('Run CellBender to remove ambient RNA?', ...
%     'py_cellbender', 'python');
% if ~ok, return; end

if ~gui.i_setpyenv([],[],FigureHandle)
    return;
end

fw=gui.myWaitbar(FigureHandle);

[output_h5] = py_cellbender(input_h5, wkdir);
% gui.myWaitbar(FigureHandle, fw, 'Processing complete. Output saved to: ', output_h5);
% pause(5);
gui.myWaitbar(FigureHandle, fw);
answer=gui.myQuestdlg(FigureHandle, sprintf('Output saved. Open the folder %s?', wkdir), '');
if strcmp(answer,'Yes'), winopen(wkdir); end

end
