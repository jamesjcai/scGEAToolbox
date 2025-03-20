function myErrordlg(parentFig, message, title)
    if nargin < 3, title = ''; end
    if nargin < 2, message = 'Message'; end
    if nargin < 1, parentFig = []; end    
    if isempty(parentFig) || ~gui.i_isuifig(parentFig)
        msgbox(message, title, 'error', 'modal');
    else
        % uialert(parentFig, message, title, 'Icon', 'error');
        uiconfirm(parentfig, message, title, ...
            'Options', {'OK'}, 'defaultoption', 'OK', ...
            'Icon', 'error');        
    end
end
