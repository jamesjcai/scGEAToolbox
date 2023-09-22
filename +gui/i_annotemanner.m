function [manuallyselect, bestonly] = i_annotemanner

manuallyselect = false;
bestonly = true;

answer = questdlg('Assign cell types to clusters automatically?', ...
    '', 'Yes, automatically', 'No, manually', ...
    'Cancel', 'Yes, automatically');
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