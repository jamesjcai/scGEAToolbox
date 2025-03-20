function h = myWarndlg(parentfig, message, title)
    if nargin < 3, title = ''; end
    if nargin < 2, message = 'Message'; end
    if nargin < 1, parentfig = []; end
    h = [];
    if isempty(parentfig) || ~gui.i_isuifig(parentfig)
        msgbox(message, title, 'warn', 'modal');
    else
        % UIFigure-based app
        uialert(parentfig, message, title, 'Icon', 'warning');        
    end
end
