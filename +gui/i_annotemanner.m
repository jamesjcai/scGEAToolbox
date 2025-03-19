function [manuallyselect, bestonly] = i_annotemanner(parentfig)

if nargin<1, parentfig = []; end

manuallyselect = false;
bestonly = true;

answer = gui.myQuestdlg(parentfig, ...
    'Assign cell types to clusters automatically?', ...
    '', {'Yes, automatically', 'No, manually', ...
    'Cancel'}, 'Yes, automatically');
if isempty(answer), return; end
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