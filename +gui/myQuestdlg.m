function answer = myQuestdlg(parentfig, message, title, options, defaultOption)
    % CUSTOMQUESTDLG Display a dialog box appropriate for the figure type.
    %
    % answer = CUSTOMQUESTDLG(parentfig, message, title, options, defaultOption)
    % displays a dialog box with the specified message and options. The behavior
    % depends on the type of parentfig:
    % - If parentfig is a traditional figure, it uses questdlg.
    % - If parentfig is a uifigure, it uses uiconfirm.
    %
    % Inputs:
    % - parentfig: Handle to the parent figure (figure or uifigure).
    % - message: Text to display in the dialog box.
    % - title: Title of the dialog box.
    % - options: Cell array of options (e.g., {'Yes', 'No', 'Cancel'}).
    % - defaultOption: Default selected option (e.g., 'Yes').
    %
    % Output:
    % - answer: The option selected by the user.


    if nargin < 4
        options = {'Yes', 'No', 'Cancel'};
        defaultOption = options{1};
    end
    if nargin < 5
        defaultOption = options{1};
    end    
    if nargin < 3, title = ''; end
    if nargin < 2, message = 'Selection'; end
    if nargin < 1, parentfig = []; end
    if isempty(parentfig) || ~gui.i_isuifig(parentfig)
        % Traditional figure-based app
        answer = questdlg(message, title, options{:}, defaultOption);
    else
        % UIFigure-based app
        answer = uiconfirm(parentfig, message, title, ...
            'Options', options, ...
            'DefaultOption', find(strcmp(options, defaultOption)), ...
            'Icon', 'question');

        if strcmp(answer, options{end})
            % if ~strcmp('Yes', gui.myQuestdlg(parentfig, ...
            %         sprintf('You selected %s. Continue?', answer)))
            %     answer = [];
            % end
            gui.myHelpdlg(parentfig, ...
                sprintf('You selected ''%s''.', answer));
        end
    end
end
