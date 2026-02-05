function callback_RunCellBender(src, ~)

[FigureHandle, ~] = gui.gui_getfigsce(src);
if ~gui.gui_showrefinfo('CellBender [PMID:37550580]', FigureHandle), return; end


if ~gui.i_setpyenv([],[],FigureHandle)
    return;
end
v = pyrun("import torch; v = torch.__version__.split('+')[0]", "v");
v = string(v);
tf = isVersionInRange(v,"1.13","2.1");
if ~tf
    % error('Unsupported PyTorch version: %s. Please use a version between 1.13 and 2.1.', v);    
    gui.myWarndlg(FigureHandle, ...
        sprintf('Unsupported PyTorch version: %s. Please use a version between 1.13 and 2.1.', v),...
        'Runtime Error', true);
    retrun;
end


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


fw=gui.myWaitbar(FigureHandle, [], false, 'Running CellBender...');

try
    [output_h5] = run.py_cellbender(input_h5, wkdir);
catch ME
    gui.myWaitbar(FigureHandle, fw, true);
    return;
end

if ~isempty(output_h5)
    gui.myWaitbar(FigureHandle, fw, [], sprintf('Processing complete. Output saved to: %s', output_h5));
    
    answer=gui.myQuestdlg(FigureHandle, sprintf('Output saved to %s. Open the folder %s?', ...
        output_h5, wkdir), '');
else
    gui.myWarndlg(FigureHandle, ...
        'CellBender terminated due to an unexpected error.',...
        'Runtime Error', true);
    answer=gui.myQuestdlg(FigureHandle, 'Open working folder to check the log output for additional details?');
end
    if strcmp(answer,'Yes'), winopen(wkdir); end    

end



function tf = isVersionInRange(ver, low, high)
% Returns true if low <= ver <= high

    v  = parseVersion(ver);
    v1 = parseVersion(low);
    v2 = parseVersion(high);

    tf = compareVersion(v, v1) >= 0 && compareVersion(v, v2) <= 0;
end


function v = parseVersion(s)
    % Convert "2.8.0" → [2 8 0]
    v = sscanf(char(s), '%d.%d.%d').';
    v(end+1:3) = 0;   % pad missing parts
end


function c = compareVersion(a,b)
    % Returns:
    %  1 if a>b
    %  0 if a=b
    % -1 if a<b

    for i = 1:3
        if a(i) > b(i)
            c = 1;
            return
        elseif a(i) < b(i)
            c = -1;
            return
        end
    end
    c = 0;
end





