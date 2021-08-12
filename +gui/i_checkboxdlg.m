function [checked] = i_checkboxdlg(options)
% determine the number of handles
if nargin<1, options={'option 1', 'option 2', 'option 3'}; end
nOptions = numel(options);
boxHeight = nOptions * 50;
boxWidth = 500;
checkBoxHeights = linspace(30, boxHeight-30, nOptions);
checkBoxHeights = flip(checkBoxHeights);
% disp('Check box heights')
% disp(checkBoxHeights)
% Create figure
h.f = figure('units','pixels','position',[400,400,boxWidth,boxHeight],...
             'toolbar','none','menu','none');
movegui(h.f, 'center');
% Create checkboxes 
for op = 1:nOptions
    h.c(op) = uicontrol('style','checkbox','units','pixels',...
                'position',[10,checkBoxHeights(op),200,15],'string',options{op});
end
% Create OK pushbutton   
h.p = uicontrol('style','pushbutton','units','pixels',...
                'position',[220,5,70,20],'string','OK',...
                'callback',@p_call);
    % Pushbutton callback
    function checked = p_call(varargin)
        vals = get(h.c,'Value');
        checked = find([vals{:}]);
        close(h.f)
        if isempty(checked)
            checked = 'none';
        end
        disp(checked)
    end
end

% https://www.mathworks.com/matlabcentral/answers/13351-dialog-with-checkboxes-in-gui
