function answer = myQuestdlg(parentfig, message, title, ...
                options, defaultOption, icontag)
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
% - answer: The option selected by the user, or '' if the dialog was
%   dismissed -- Escape, the window's close button, or a Cancel button this
%   function added itself. That matches QUESTDLG, whose '' on dismissal is
%   what the 307 call sites in this toolbox test for with
%   `if isempty(answer), return; end`.
%
%   A caller that lists 'Cancel' among its own OPTIONS is branching on the
%   label and gets it back verbatim; 50 of the 307 do, GUI.I_TRANSFORMX
%   among them.

if nargin < 6 || isempty(icontag)
    icontag = 'question'; % warning
end
if nargin < 4 || isempty(options)
    options = {'Yes', 'No', 'Cancel'};
    defaultOption = options{1};
end

if nargin < 5 || isempty(defaultOption)   % put this here afte argin<4
    defaultOption = options{1};
end

if nargin < 3, title = ''; end
if nargin < 2, message = 'Selection'; end
if nargin < 1, parentfig = []; end

if ~isempty(parentfig) && pkg.i_isvalid(parentfig) && parentfig.Visible == "on"
    figure(parentfig);
    % The restore on the way out goes through GUI.I_RAISEFIG rather than a
    % bare FIGURE(): the handle can be gone by the time the cleanup runs --
    % a callback that closes the main window while the dialog is up -- and a
    % destructor that throws surfaces as a warning the user cannot act on.
    % It also raises a uifigure properly, which FIGURE() does not.
    cleanupObj = onCleanup(@() gui.i_raisefig(parentfig));
end

if isempty(parentfig) || ~gui.i_isuifig(parentfig)
    % Traditional figure-based app
    answer = questdlg(message, title, options{:}, defaultOption);
else
    % UIFigure-based app
    %
    % UICONFIRM has no equivalent of QUESTDLG's empty return: it needs a
    % CancelOption, and hands back that option's LABEL when the user
    % presses Escape, closes the window or clicks it. So this branch used
    % to return the literal 'Cancel' to callers that were watching for ''
    % -- and scgeatool is App Designer-based, so this is the branch the
    % main GUI takes. Every `if isempty(answer), return; end` therefore
    % failed to abort. gui.callback_Harmony, callback_Harmonypy and
    % callback_HarmonyR each went on to run a backend the user had just
    % cancelled, and callback_Harmonypy's dialog defaults to the Python
    % path, so pressing Escape started a Python run.
    %
    % Only the Cancel THIS FUNCTION adds is normalised. A caller that
    % lists 'Cancel' among its own options is branching on the label and
    % must keep receiving it: gui.i_transformx does exactly that, and
    % raises 'Wrong option' on anything it does not recognise.
    cancelWasOurs = ~any(strcmp(options, 'Cancel'));
    if cancelWasOurs
        options{end+1} = 'Cancel';
    end

    answer = uiconfirm(parentfig, message, title, ...
        'Options', options, ...
        'DefaultOption', find(strcmp(options, defaultOption)), ...
        'Icon', icontag, 'CancelOption', length(options));

    if cancelWasOurs && strcmp(answer, 'Cancel')
        answer = '';   % dismissed, reported as questdlg reports it
    end

    % if strcmp(answer, options{end})
    %     % if ~strcmp('Yes', gui.myQuestdlg(parentfig, ...
    %     %         sprintf('You selected %s. Continue?', answer)))
    %     %     answer = [];
    %     % end
    %     gui.myHelpdlg(parentfig, ...
    %         sprintf('You selected ''%s''.', answer));
    % end
end
end
