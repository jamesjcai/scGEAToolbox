function answer = myQuestdlg(parentFig, message, title, options, defaultOption)
    % CUSTOMQUESTDLG Display a dialog box appropriate for the figure type.
    %
    % answer = CUSTOMQUESTDLG(parentFig, message, title, options, defaultOption)
    % displays a dialog box with the specified message and options. The behavior
    % depends on the type of parentFig:
    % - If parentFig is a traditional figure, it uses questdlg.
    % - If parentFig is a uifigure, it uses uiconfirm.
    %
    % Inputs:
    % - parentFig: Handle to the parent figure (figure or uifigure).
    % - message: Text to display in the dialog box.
    % - title: Title of the dialog box.
    % - options: Cell array of options (e.g., {'Yes', 'No', 'Cancel'}).
    % - defaultOption: Default selected option (e.g., 'Yes').
    %
    % Output:
    % - answer: The option selected by the user.

%    if nargin < 5
%        defaultOption = options{1};
%    end
    if nargin < 4
        options = {'Yes', 'No', 'Cancel'};
        defaultOption = options{1};
    end
    if nargin < 3, title = ''; end
    if nargin < 2, message = 'Selection'; end
    if nargin < 1, parentFig = []; end
    if isempty(parentFig) || ~gui.i_isuifig(parentFig)
        % Traditional figure-based app
        answer = questdlg(message, title, options{:}, defaultOption);
    else
        % UIFigure-based app
        answer = uiconfirm(parentFig, message, title, ...
            'Options', options, ...
            'DefaultOption', find(strcmp(options, defaultOption)), ...
            'Icon', 'question');
    end
end
