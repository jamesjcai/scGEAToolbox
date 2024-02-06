% function yourDlg(Str, List)
% dialog('Position', [100, 100, 400, 600]);
% uicontrol('Style', 'ListBox', 'String', List, ...
%           'Position', [10, 530, 380, 30]);
% uicontrol('Style', 'Edit', 'String', Str, ...
%           'Position', [10, 10, 380, 500]);
% end

formats = struct('type', {}, 'style', {}, 'items', {}, ...
  'format', {}, 'limits', {}, 'size', {});
formats(1,1).type   = 'edit';
formats(1,1).format = 'integer';
formats(1,1).limits = [1 50];

formats(2,1).type   = 'list';
formats(2,1).style  = 'popupmenu';
formats(2,1).items  = {'jet', 'hsv', 'hot', 'cool', 'spring', 'summer', ...
  'autumn', 'winter', 'gray', 'bone', 'copper', 'pink', 'lines'};
defaultanswer = {20, 2};

[answer, canceled] = inputsdlg('aaa', 'bbb', formats, defaultanswer);