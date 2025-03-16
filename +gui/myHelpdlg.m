function [h] = myHelpdlg(parentfig, message, title)
    if nargin < 3, title = ''; end
    if nargin < 2, message = 'Message'; end
    if nargin < 1, parentfig = []; end
    if nargout > 0, h = []; end
    
    if isempty(parentfig) || ~gui.i_isuifig(parentfig)
        % Traditional figure-based app
        % h = helpdlg(message, title);
        msgbox(message, title, 'help', 'modal');
    else
        % UIFigure-based app
        % disp('aaa UIFigure-based app')
        % uialert(parentfig, message, title, 'Icon', 'info');
        uiconfirm(parentfigx, message, title, ...
            'Options', {'OK'}, 'defaultoption', 'OK', ...
            'Icon', 'info');
    end
end
