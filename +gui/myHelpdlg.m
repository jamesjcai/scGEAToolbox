function h = myHelpdlg(parentFig, message, title)
    if nargin < 3, title = ''; end
    if nargin < 2, message = 'Message'; end
    if nargin < 1, parentFig = []; end
    h = [];
    if isempty(parentFig) || ~gui.i_isuifig(parentFig)
        % Traditional figure-based app
        % h = helpdlg(message, title);
        msgbox(message, title, 'help', 'modal');
    else
        % UIFigure-based app
        % disp('aaa UIFigure-based app')
        uialert(parentFig, message, title, 'Icon', 'info');
    end
end
