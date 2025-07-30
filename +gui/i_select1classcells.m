function [ptsSelected] = i_select1classcells(sce, askunselect, parentfig)
if nargin < 2, askunselect = true; end
if nargin < 3, parentfig = []; end
ptsSelected = [];


        [thisc, clabel] = gui.i_select1state(sce, ...
            false, false, true, false, parentfig);

%[thisc, clabel] = gui.i_select1class(sce,[],[],[],parentfig);
if isempty(thisc), return; end

answer2 = gui.myQuestdlg(parentfig, sprintf('How to sort members of ''%s''?',clabel), '', ...
    {'Alphabetic', 'Size (Descending Order)'}, 'Alphabetic');
switch answer2
    case 'Alphabetic'
        [~, cLi] = findgroups(string(thisc));
        [cLisorted] = natsort(string(cLi));
    case 'Size (Descending Order)'
        [cLisorted]=pkg.e_sortcatbysize(thisc);
    otherwise
        return;
end

if gui.i_isuifig(parentfig)
    [indxx, tfx] = gui.myListdlg(parentfig, cLisorted, 'Select groups');
else
    [indxx, tfx] = listdlg('PromptString', {'Select groups'}, ...
        'SelectionMode', 'multiple', ...
        'ListString', cLisorted, ...
        'ListSize', [220, 300]);
end

if tfx == 1
    ptsSelected = ismember(string(thisc), cLisorted(indxx));
    %ptsSelected=ismember(ci,indxx);
    if askunselect
        answer = gui.myQuestdlg(parentfig, 'Select or unselect?', '', {'Select', 'Unselect', ...
            'Cancel'}, 'Select');
        if strcmp(answer, 'Select')
        elseif strcmp(answer, 'Unselect')
            ptsSelected = ~ptsSelected;
        else
            ptsSelected = [];
        end
    end
end
end
