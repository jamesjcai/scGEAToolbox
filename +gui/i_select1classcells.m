function [ptsSelected] = i_select1classcells(sce, askunselect)
if nargin < 2, askunselect = true; end
ptsSelected = [];
thisc = gui.i_select1class(sce);
if isempty(thisc), return; end

[~, cLi] = grp2idx(thisc);

[cLisorted] = natsort(string(cLi));
[indxx, tfx] = listdlg('PromptString', {'Select groups'}, ...
    'SelectionMode', 'multiple', 'ListString', cLisorted);
if tfx == 1
    ptsSelected = ismember(string(thisc), cLisorted(indxx));
    %ptsSelected=ismember(ci,indxx);
    if askunselect
        answer = questdlg('Select or unselect?', '', 'Select', 'Unselect', ...
            'Cancel', 'Select');
        if strcmp(answer, 'Select')
        elseif strcmp(answer, 'Unselect')
            ptsSelected = ~ptsSelected;
        else
            ptsSelected = [];
        end
    end
end
end
