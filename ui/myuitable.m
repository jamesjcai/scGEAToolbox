% https://stackoverflow.com/questions/62456301/how-can-i-make-a-matlab-gui-that-takes-many-input-values-and-can-be-displayed-wi
prompt = {'Enter a value:',... %answer 1
    'Enter a value:',...       %answer 2
    'Enter a value:',...       %answer 3
    'Enter a value:',...       %answer 4
    'Enter a value:',...       %answer 5
    'Enter a value:',...       %answer 6
    'Enter a value:',...       %answer 7
    'Enter a value:',...       %answer 8
    'Enter a value:',...       %answer 9
    'Enter a value:',...       %answer 10
    'Enter a value:',...       %answer 11
    'Enter a value:',...       %answer 12
    'Enter a value:',...       %answer 13
    'Enter a value:',...       %answer 14
    'Enter a value:',...       %answer 15
    'Enter a value:',...       %answer 16
    'Enter a value:',...       %answer 17
    'Enter a value:',...       %answer 18
    };
title = 'Specifications';            
dims = [1 35];              % input field specifications
definput = {'0','0','1','0.5','2.0','15','0.3','1','1','2','1.0','20','0','3000','2^7','2','0','Y'};   % default values

% Create figure & uitable
fig = uifigure('Name',title,'Position',[100 100 500 500],...
    'CloseRequestFcn',@(fig,event) Cancel_Pushed(fig)); % Create figure, set window size and close request
tableData = {prompt{:};definput{:}}'; % Data for table
table = uitable(fig); % Create table
table.Position = [50 80 400 400]; % Set table size
table.Data = tableData; % Send data to table
table.ColumnEditable = [false,true]; % Allow input on second column

% Create buttons w/ callbacks
OK = uibutton(fig,'Position',[50 30 100 22],'Text','OK','ButtonPushedFcn', @(OK,event) OK_Pushed(fig,table));
Cancel = uibutton(fig,'Position',[170 30 100 22],'Text','Cancel','ButtonPushedFcn', @(Cancel,event) Cancel_Pushed(fig));

dataOut = []; % Declare global variable to retrieve data from callbacks
global dataOut; 

uiwait(fig); % Wait until figure is closed


function OK_Pushed(fig,table)
    global dataOut
    tableOut = table.Data; % Get data
    dataOut = tableOut(:,2); % Save table data
    delete(fig); % Close figure

end
function Cancel_Pushed(fig)
    global dataOut
    dataOut = {0}; % Return 0 if canceled
    delete(fig); % Close figure
end