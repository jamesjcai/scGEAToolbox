function [h] = myErrordlg(parentfig, message, title, modal)

    if nargin < 4, modal = true; end
    if nargin < 3, title = ''; end
    if nargin < 2, message = 'Message'; end
    if nargin < 1, parentfig = []; end
    gui.mydlg('error', parentfig, message, title, modal);

    %{
    if nargin < 4, modal = true; end
    if nargin < 3, title = ''; end
    if nargin < 2, message = 'Message'; end
    if nargin < 1, parentfig = []; end
    if nargout > 0, h = []; end
    if isempty(parentfig) || ~gui.i_isuifig(parentfig)
        if modal
            msgbox(message, title, 'error', 'modal');
        else
            msgbox(message, title, 'error');
        end            
    else
        if modal
            h = uiconfirm(parentfig, message, title, ...
                'Options', {'OK'}, 'defaultoption', 'OK', ...
                'Icon', 'error');
        else
            uialert(parentfig, message, title, 'Icon', 'error');
            %uiconfirm(parentfig, message, title, ...
            %    'Options', {'OK'}, 'defaultoption', 'OK', ...
            %    'Icon', 'error');        
        end
    end
end
    %}