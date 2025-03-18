function [newclabel] = i_renamec(clabel, sce, newpickclabel, parentfig)
if nargin < 4, parentfig = []; end
if nargin < 3, newpickclabel = ''; end

if ~isempty(newpickclabel)
    newclabel = newpickclabel;
else
    newclabel = matlab.lang.makeValidName('Workspace Variable');
end

if strcmp(clabel, 'Workspace Variable...')
    answerxx = gui.myQuestdlg(parentfig, 'Rename Workspace Variable?','');
    if strcmp(answerxx, 'Yes')


if gui.i_isuifig(parentfig)
    newclabelx = gui.myInputdlg({'New name'}, 'Rename', {newpickclabel}, parentfig);
else
    newclabelx = inputdlg('New name', 'Rename', [1, 50], {newpickclabel});
end

        if ~isempty(newclabelx) && ~isempty(newclabelx{1})
            newclabel = matlab.lang.makeValidName(newclabelx);
        end
    end
    currentclabel = sce.list_cell_attributes(1:2:end);
    currentclabelx = [currentclabel, newclabel];
    currentclabelx1 = matlab.lang.makeUniqueStrings(currentclabelx);
    newclabel = currentclabelx1(end);
    if ~isequal(currentclabelx1, currentclabelx)
        gui.myHelpdlg(parentfig, sprintf(['Name exisiting. New name ' ...
            'is changed to %s.'], ...
            newclabel{1}));
    end
end
end
