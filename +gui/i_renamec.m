function [newclabel] = i_renamec(clabel, sce, newpickclabel)

if nargin < 3, newpickclabel = ''; end

if ~isempty(newpickclabel)
    newclabel = newpickclabel;
else
    newclabel = matlab.lang.makeValidName('Workspace Variable');
end

if strcmp(clabel, 'Workspace Variable...')
    answerxx = questdlg('Rename Workspace Variable?');
    if strcmp(answerxx, 'Yes')
        newclabelx = inputdlg('New name', ...
            'Rename', [1, 50], {newpickclabel});
        if ~isempty(newclabelx) && ~isempty(newclabelx{1})
            newclabel = matlab.lang.makeValidName(newclabelx);
        end
    end
    currentclabel = sce.list_cell_attributes(1:2:end);
    currentclabelx = [currentclabel, newclabel];
    currentclabelx1 = matlab.lang.makeUniqueStrings(currentclabelx);
    newclabel = currentclabelx1(end);
    if ~isequal(currentclabelx1, currentclabelx)
        waitfor(helpdlg(sprintf('Name exisiting. New name is changed to %s.', ...
            newclabel{1}), ''));
    end
end
end
