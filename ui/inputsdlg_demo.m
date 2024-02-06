%INPUTSDLG DEMO (Enhanced input dialog box with multiple data types)

% Written by: Takeshi Ikuma
% Last Updated: May 5 2010
%
% Updated for additional functions F. Hatz 2013

clear; close all;

Title = 'INPUTSDLG Demo Dialog';

%%%% SETTING DIALOG OPTIONS
% Options.WindowStyle = 'modal';
Options.Resize = 'on';
Options.Interpreter = 'tex';
Options.CancelButton = 'on';
Options.ApplyButton = 'on';
Options.ButtonNames = {'Continue','Cancel'}; %<- default names, included here just for illustration
Option.Dim = 4; % Horizontal dimension in fields

Prompt = {};
Formats = {};
DefAns = struct([]);

Prompt(1,:) = {['This demo illustrates every type of control that can be placed by INPUTSDLG function ' ... 
   'and demonstrates how Formats input can be used to layout these controls.'],[],[]};
Formats(1,1).type = 'text';
Formats(1,1).size = [-1 0];
Formats(1,1).span = [1 2]; % item is 1 field x 4 fields

Prompt(2,:) = {'Bidder''s Name', 'Name',[]};
Formats(2,1).type = 'edit';
Formats(2,1).format = 'text';
Formats(2,1).size = 200; % automatically assign the height
DefAns(1).Name = 'John Smith';

Prompt(3,:) = {'Bidder''s SSN', 'SSN','(no space or hyphen)'};
Formats(2,2).type = 'edit';
Formats(2,2).format = 'integer';
Formats(2,2).limits = [0 999999999]; % 9-digits (positive #)
Formats(2,2).size = 80;
Formats(2,2).unitsloc = 'bottomleft';
DefAns.SSN = 123456789;

Prompt(4,:) = {'Bidding Price: $', 'Price',[]};
Formats(2,3).type = 'edit';
Formats(2,3).format = 'float';
Formats(2,3).margin = [0 0];
Formats(2,3).limits = [0 inf]; % non-negative decimal number
Formats(2,3).callback = @(hObj,~,~,~)set(hObj,'String',num2str(round(str2double(get(hObj,'String'))*100)/100));
DefAns.Price = 99.99;

Prompt(5,:) = {'Enable enhanced mode' 'EnableEnhancedMode',[]};
Formats(2,4).type = 'check';
DefAns.EnableEnhancedMode = true;

Prompt(end+1,:) = {'Bidder''s Bio File','BioFile',[]};
Formats(3,1).type = 'edit';
Formats(3,1).format = 'file';
Formats(3,1).items = {'*.bio','Biography File (*.bio)';'*.*','All Files'};
Formats(3,1).limits = [0 1]; % single file get
Formats(3,1).size = [-1 0];
Formats(3,1).span = [1 3];  % item is 1 field x 3 fields
d = dir;
files = strcat([pwd filesep],{d(~[d.isdir]).name});
DefAns.BioFile = files{1};

Prompt(end+1,:) = {'Action','Action',[]};
Formats(3,4).type = 'list';
Formats(3,4).style = 'togglebutton';
Formats(3,4).format = 'text';
Formats(3,4).items = {'Bid';'Decline';'Pass'};
Formats(3,4).span = [3 1];  % item is 3 fields x 1 field
DefAns.Action = 'Decline'; % = 'Decline'

Prompt(end+1,:) = {'Bidder''s Data Folder','DataFolder',[]};
Formats(4,1).type = 'edit';
Formats(4,1).format = 'dir';
Formats(4,1).size = [-1 0];
Formats(4,1).span = [1 3];  % item is 1 field x 3 fields
DefAns.DataFolder = pwd;

Prompt(end+1,:) = {'Save Bidding History To','SaveFile',[]};
Formats(5,1).type = 'edit';
Formats(5,1).format = 'file';
Formats(5,1).limits = [1 0]; % use uiputfile
Formats(5,1).size = [-1 0];
Formats(5,1).span = [1 3];  % item is 1 field x 3 fields
DefAns.SaveFile = files{2};

Prompt(end+1,:) = {'Select Item Files','ItemFiles',[]};
Formats(6,1).type = 'edit';
Formats(6,1).format = 'file';
Formats(6,1).limits = [0 5]; % multi-select files
Formats(6,1).size = [-1 -1];
Formats(6,1).items = {'*.itm','Auction Item File';'*.*','All Files'};
Formats(6,1).span = [1 3];  % item is 1 field x 3 fields
DefAns.ItemFiles = files(3:end);

Prompt(end+1,:) = {'Choose Currency','MoneyUnit',[]};
Formats(7,1).type = 'list';
Formats(7,1).format = 'text';
Formats(7,1).style = 'radiobutton';
Formats(7,1).items = {'U.S. Dollar' 'Euro';'Japanese Yen' ''};
% Formats(7,1).span = [2 1];  % item is 2 field x 1 fields
DefAns.MoneyUnit = 'Japanese Yen';%3; % yen

Prompt(end+1,:) = {'Item Table','Table',[]};
Formats(7,2).type = 'table';
Formats(7,2).format = {'char', {'left','right'}, 'numeric' 'logical'}; % table (= table in main dialog) / window (= table in separate dialog)
Formats(7,2).items = {'Row' 'Left/Right' 'Value' 'Done'};
Formats(7,2).size = [373 73];
Formats(7,2).span = [1 3];  % item is 2 field x 1 fields
DefAns.Table = {'Row1' 'left'  10 false
                'Row2' 'right' 12 true};

Prompt(end+1,:) = {'Bidding Rate','BidRate',[]};
Formats(8,1).type = 'range';
DefAns.BidRate = 0.75;

Prompt(end+1,:) = {'Press Me!','',''};
Formats(8,2).type = 'button';
Formats(8,2).callback = @(~,~,handles,k)msgbox(sprintf('You just pressed %s button',get(handles(k),'String')),'modal');

Prompt(end+1,:) = {'Item List','List',[]};
Formats(8,3).type = 'list';
Formats(8,3).style = 'popupmenu';
Formats(8,3).items = {'Black','White','Red','Blue','Green','Yellow','Orange'};

Prompt(end+1,:) = {'Item Color','Color',[]};
Formats(8,4).type = 'color';
DefAns.Color = [1 0 0];

Prompt(end+1,:) = {'Item Vector','Vector',[]};
Formats(9,1).type = 'edit';
Formats(9,1).format = 'vector';
Formats(9,1).limits = [0 5]; % 5 parameters max

Prompt(end+1,:) = {'Date (dd-mmm-yyyy)', 'Date',[]};
Formats(9,2).type = 'edit';
Formats(9,2).format = 'date';
Formats(9,2).limits = 'dd-mmm-yyyy'; % with time: 'dd-mmm-yyyy HH:MM:SS'
Formats(9,2).size = 100;
DefAns.Date = [2013 1 2 0 0 0];

Prompt(end+1,:) = {'Memo','Memo',[]};
Formats(10,1).type = 'edit';
Formats(10,1).format = 'text';
Formats(10,1).limits = [0 9]; % default: show 20 lines
Formats(10,1).size = [-1 0];
Formats(10,1).span = [1 4];  % item is 4 fields x 3 fields

Prompt(end+1,:) = {'Auction Sites:','Site',[]};
Formats(11,1).type = 'list';
Formats(11,1).style = 'listbox';
Formats(11,1).format = 'text'; % Answer will give value shown in items, disable to get integer
Formats(11,1).items = {'www.auction1.com';'www.auction2.com';'www.bidme.com';'www.bestvalu.com'};
Formats(11,1).limits = [0 4]; % multi-select
Formats(11,1).size = [140 80];
% Formats(11,1).span = [3 1];  % item is 2 fields x 1 fields
DefAns.Site = {'www.auction1.com','www.bidme.com'};

Prompt(end+1,:) = {'X', 'X','(-4, 12)'};
Formats(11,2).type = 'edit';
Formats(11,2).format = 'float';
Formats(11,2).limits = [-4 12];
Formats(11,2).size = 80;
Formats(11,2).callback = @(~,~,h,k)set(h(k+2),'String',num2str(str2double(get(h(k),'String'))*str2double(get(h(k+1),'String'))));

Prompt(end+1,:) = {'Y', 'Y','(non-negative)'};
Formats(11,3).type = 'edit';
Formats(11,3).format = 'integer';
Formats(11,3).limits = [0 100];
Formats(11,3).size = 80;
Formats(11,3).callback = @(~,~,h,k)set(h(k+1),'String',num2str(str2double(get(h(k-1),'String'))*str2double(get(h(k),'String'))));

Prompt(end+1,:) = {'X*Y','XY',[]};
Formats(11,4).type = 'edit';
Formats(11,4).format = 'float';
Formats(11,4).enable = 'inactive';
Formats(11,4).size = [100 25];

DefAns.X = pi;
DefAns.Y = 4;
DefAns.XY = DefAns.X*DefAns.Y;

[Answer,Cancelled] = inputsdlg(Prompt,Title,Formats,DefAns,Options)
