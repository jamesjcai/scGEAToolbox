function newclable=i_renamec(clable,sce,newpickclable)

if nargin<3, newpickclable=''; end

if ~isempty(newpickclable)
    newclable=newpickclable;
else
    newclable=matlab.lang.makeValidName('Customized C');
end

if strcmp(clable,'Customized C...')
    answerxx = questdlg('Rename Customized C?');
    if strcmp(answerxx, 'Yes')
        newclablex = inputdlg('New name',...
            'Rename', [1 50], {newpickclable});
        if ~isempty(newclablex) && ~isempty(newclablex{1})
            newclable=matlab.lang.makeValidName(newclablex);
        end
    end
    currentclable=sce.list_cell_attributes(1:2:end);
    currentclablex=[currentclable, newclable];
    currentclablex1 = matlab.lang.makeUniqueStrings(currentclablex);
    newclable=currentclablex1(end);
    if ~isequal(currentclablex1,currentclablex)
        waitfor(helpdlg(sprintf('Name exisiting. New name is changed to %s.',...
            newclable{1}),''));
    end
end
end

