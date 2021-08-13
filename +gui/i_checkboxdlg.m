function [checked] = i_checkboxdlg(options)
% determine the number of handles
if nargin<1, options={'Number of cells > 500',...
                      'Number of genes < 10,000',...
                      'option 3'};
end
nOptions = numel(options);
boxHeight = nOptions * 33;
boxWidth = 310;
checkBoxHeights = linspace(30, boxHeight-30, nOptions);
checkBoxHeights = flip(checkBoxHeights);
% disp('Check box heights')
% disp(checkBoxHeights)
% Create figure
% h.f = figure('units','pixels','position',[400,400,boxWidth,boxHeight],...
%              'toolbar','none','menu','none','Name','Requirement Checklist',...
%              'NumberTitle','off');
% movegui(h.f, 'center');

h.f = dialog('Position',[400,400,boxWidth,boxHeight],...
    'Name','Requirement Checklist','Visible',false);
movegui(h.f, 'center');
set(h.f,'visible',true);

% Create checkboxes 
for op = 1:nOptions
    h.c(op) = uicontrol('Parent',h.f,'style','checkbox','units','pixels',...
                'position',[10,checkBoxHeights(op),200,15],...
                'string',options{op},'Value',1);
end
checked=[];
% Create OK pushbutton   
h.p = uicontrol('Parent',h.f,'style','pushbutton','units','pixels',...
                'position',[130,10,70,20],'string','OK',...
                'callback',@p_call);
        %vals = get(h.c,'Value');
        %checked = find([vals{:}]);
            
    % Pushbutton callback
    function checked = p_call(varargin)
        vals = get(h.c,'Value');
        checked = find([vals{:}]);
        close(h.f)
        %if isempty(checked)
        %    checked = 'none';
        %end
        %disp(checked)
    end
end

% https://www.mathworks.com/matlabcentral/answers/13351-dialog-with-checkboxes-in-gui
