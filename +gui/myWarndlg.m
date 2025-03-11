function h = myWarndlg(parentFig, message, title)
    if nargin < 3, title = ''; end
    if nargin < 2, message = 'Message'; end
    if nargin < 1, parentFig = []; end
    h = [];
    if isempty(parentFig) || ~gui.i_isuifig(parentFig)
        msgbox(message, title, 'warn', 'modal');
    else
        % UIFigure-based app
        uialert(parentFig, message, title, 'Icon', 'warning');        
    end
end
