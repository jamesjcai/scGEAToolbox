function myErrordlg(parentfig, message, title)
    if nargin < 3, title = ''; end
    if nargin < 2, message = 'Message'; end
    if nargin < 1, parentfig = []; end    
    if isempty(parentfig) || ~gui.i_isuifig(parentfig)
        msgbox(message, title, 'error', 'modal');
    else
        % uialert(parentfig, message, title, 'Icon', 'error');
        uiconfirm(parentfig, message, title, ...
            'Options', {'OK'}, 'defaultoption', 'OK', ...
            'Icon', 'error');        
    end
end
