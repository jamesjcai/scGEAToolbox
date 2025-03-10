function h = myWarndlg(parentFig, message, title)
    if nargin < 3, title = ''; end
    if nargin < 2, message = 'Message'; end
    if nargin < 1, parentFig = []; end
    if isempty(parentFig) || ~gui.i_isuifig(parentFig)
        % Traditional figure-based app
        h = warndlg(message, title);
    else
        % UIFigure-based app
        uialert(parentFig, message, title, 'Icon', 'warning');
        h = [];
    end
end
