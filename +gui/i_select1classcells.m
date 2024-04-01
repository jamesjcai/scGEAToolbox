function [ptsSelected] = i_select1classcells(sce, askunselect)
if nargin < 2, askunselect = true; end
ptsSelected = [];
[thisc, clabel] = gui.i_select1class(sce);
if isempty(thisc), return; end



answer2 = questdlg(sprintf('How to sort members of ''%s''?',clabel), '', ...
    'Alphabetic', 'Size (Descending Order)', 'Alphabetic');
switch answer2
    case 'Alphabetic'
        [~, cLi] = grp2idx(thisc);
        [cLisorted] = natsort(string(cLi));
    case 'Size (Descending Order)'
        [cLisorted]=pkg.e_sortcatbysize(thisc);
    otherwise
        return;
end

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
