function [idx] = i_selmultidlg(genelist, predefinedlist)
%gui.i_selmultidlg
% This function is called by gui.i_selectngenes

    parentfig = gcf;
    p = parentfig.Position;
    cx = [p(1)+p(3)/2 p(2)+p(4)/2];


idx = [];
if nargin < 2, predefinedlist = []; end
txt_cell_array = {'line1'; 'line2'; 'line3'; 'line4'; 'line5'};
if nargin < 1, genelist = txt_cell_array; end
if ~iscell(genelist)
    genelist = cellstr(genelist);
end
if ~isempty(predefinedlist)
    if ~iscell(predefinedlist)
        %predefinedlist=cellstr(predefinedlist);
    end
    %predefinedlist=genelist(matches(genelist,predefinedlist,'IgnoreCase',true));
end
if ~isempty(predefinedlist)
    inlist = setxor(genelist, predefinedlist, 'stable');
else
    inlist = genelist;
end



f = figure('Visible', 'off', ...
    'WindowStyle', 'modal', 'NumberTitle', 'off');
% ax = axes(f);
% ax.Units = 'pixels';
% ax.Position = [75 75 325 280]

    px = f.Position;
    px_new = [cx(1)-px(3)/2 cx(2)-px(4)/2];
    movegui(f, px_new);

% movegui(f, 'center');
uicontrol('style', 'pushbutton', 'Position', [220, 180, 100, 30], ...
    'String', '>', 'Callback', {@plotButtonPushed, genelist});
uicontrol('style', 'pushbutton', 'Position', [220, 215, 100, 30], ...
    'String', '<', 'Callback', {@plotButtonPushed2, genelist});
uicontrol('style', 'pushbutton', 'Position', [220, 125, 100, 30], ...
    'String', 'Done', 'Callback', {@doneButtonPushed, genelist});


h_list1 = uicontrol('style', 'list', 'max', length(genelist), ...
    'min', 1, 'Position', [20, 20, 170, 360], ...
    'string', inlist);

h_list2 = uicontrol('style', 'list', 'max', length(genelist), ...
    'min', 1, 'Position', [360, 20, 170, 360], ...
    'string', predefinedlist);
set(f, 'Visible', 'on')
drawnow();
uiwait(); % add this to the end

    function plotButtonPushed(~, ~, genelist)
        if ~isempty(h_list1.String)
            h_list2.String = ...
                unique([h_list2.String; ...
                h_list1.String(h_list1.Value)], 'stable');
            h_list1.String = setxor(genelist, h_list2.String, 'stable');
            set(h_list1, 'Value', 1);
        end
end

   function plotButtonPushed2(~, ~, genelist)
        % https://www.mathworks.com/matlabcentral/answers/92064-why-do-i-receive-a-warning-when-i-repopulate-my-listbox-uicontrol-in-matlab
        % set(h_list2,'Value',1);
        if ~isempty(h_list2.String)
            h_list2.String(h_list2.Value) = [];
            set(h_list2, 'Value', 1);
        end
        h_list1.String = setxor(genelist, h_list2.String, 'stable');
    end

    function doneButtonPushed(~, ~, genelist)
        [~, idx] = ismember(h_list2.String, genelist);
        uiresume();
        closereq();
    end
end

