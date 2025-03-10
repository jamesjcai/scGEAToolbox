function [manuallyselect, bestonly] = i_annotemanner(FigureHandle)

if nargin<1, FigureHandle = []; end

manuallyselect = false;
bestonly = true;

answer = gui.myQuestdlg(FigureHandle, 'Assign cell types to clusters automatically?', ...
    '', {'Yes, automatically', 'No, manually', ...
    'Cancel'}, 'Yes, automatically');
switch answer
    case 'Yes, automatically'
        manuallyselect = false;
        bestonly = true;
    case 'No, manually'
        manuallyselect = true;
        bestonly = false;
    otherwise
        return;
end
end