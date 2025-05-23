function mydlg(tag, parentfig, message, title, modal)
    if nargin < 5, modal = true; end
    if nargin < 4, title = ''; end
    if nargin < 3, message = 'Message'; end
    if nargin < 2, parentfig = []; end
    if nargin < 1, tag = 'info'; end    % info error warning
    if isempty(parentfig) || ~gui.i_isuifig(parentfig)
        if modal
            waitfor(msgbox(message, title, tag, 'modal'));  % This will pause execution until the msgbox is closed 
        else
            msgbox(message, title, tag, 'modal');
        end            
    else
        if modal
            answer = uiconfirm(parentfig, message, title, ...
                'Options', {'OK'}, 'defaultoption', 'OK', ...
                'Icon', tag);
        else
            uialert(parentfig, message, title, 'Icon', tag);
        end
    end
end
