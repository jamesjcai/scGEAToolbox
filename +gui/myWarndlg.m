function myWarndlg(parentfig, message, title, modal)

    if nargin < 4, modal = true; end
    if nargin < 3, title = ''; end
    if nargin < 2, message = 'Message'; end
    if nargin < 1, parentfig = []; end
    if gui.i_isuifig(parentfig)
        tag = 'warning'; 
    else
        tag = 'warn'; 
    end
    gui.mydlg(tag, parentfig, message, title, modal);

%{
    if nargin < 4, modal = true; end
    if nargin < 3, title = ''; end
    if nargin < 2, message = 'Message'; end
    if nargin < 1, parentfig = []; end
    h = [];
    if isempty(parentfig) || ~gui.i_isuifig(parentfig)
        if modal
            msgbox(message, title, 'warn', 'modal');
        else
            msgbox(message, title, 'warn');
        end        
    else
        % UIFigure-based app
        if modal
            h = uiconfirm(parentfig, message, title, ...
                'Options', {'OK'}, 'defaultoption', 'OK', ...
                'Icon', 'warning');            
        else
            uialert(parentfig, message, title, 'Icon', 'warning');
        end
    end
end
%}