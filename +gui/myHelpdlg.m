function myHelpdlg(parentfig, message, title, modal)

    if nargin < 4, modal = true; end
    if nargin < 3, title = ''; end
    if nargin < 2, message = 'Message'; end
    if nargin < 1, parentfig = []; end

    if gui.i_isuifig(parentfig)
        tag = 'info'; 
    else
        tag = 'help'; 
    end

    gui.mydlg(tag, parentfig, message, title, modal);
    
    %{
    if isempty(parentfig) || ~gui.i_isuifig(parentfig)
        % Traditional figure-based app
        % h = helpdlg(message, title);
        if modal
            msgbox(message, title, 'help', 'modal');
        else
            msgbox(message, title, 'help');
        end
    else
        % UIFigure-based app
        % disp('aaa UIFigure-based app')
        % uialert(parentfig, message, title, 'Icon', 'info');
        if modal
            h = uiconfirm(parentfig, message, title, ...
                'Options', {'OK'}, 'defaultoption', 'OK', ...
                'Icon', 'info');
        else
            % uiconfirm(parentfig, message, title, ...
            %     'Options', {'OK'}, 'defaultoption', 'OK', ...
            %     'Icon', 'info');
            uialert(parentfig, message, title, 'Icon', 'info');
        end
    end
end
    %}