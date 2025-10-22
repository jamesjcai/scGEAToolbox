function [pval, tf] = sc_telemetrynotice
% setpref('scgeatoolbox', 'useronboardingtoolbar', true);

group = "scgeatoolbox";
pref =  "sharediagnosticsusage";        
title = "Help Improve SCGEATOOL";
quest = ["Share optional diagnostic data to help us improve SCGEATOOL. ", ...
         "This includes anonymous usage statistics and error reports.",...
         "No personal information is collected."];
pbtns = ["Yes","No"];

[pval,tf] = uigetpref(group,pref,title,quest,pbtns,...
    "CheckboxState",1,"DefaultButton","No");

% "ExtraOptions","Cancel");
% p = uisetpref('clearall')
% setpref(group,pref,'ask')

%{
fig = uifigure; % ("Position",[100 100 300 300]);
g = uigridlayout(fig);
%g.RowHeight =   {'1x',22};
%g.ColumnWidth = {250,100};

lbl = uilabel(g,'text','Share optional diagnostic data to help improve SCGEATOOL.');
%lbl.Layout.Row = 2;
%lbl.Layout.Column = 1;

% lbl.Position = [20 200 260 50]; % Set label position
s = uiswitch(g,"Value","On");
%s.Layout.Row = 2;
%s.Layout.Column = 2;

cbx = uicheckbox(g,"Text","Do not show this again","Value",true);

% s.Items=["",""];
b = uibutton(g,"Text","Continue");
%}


%{
function createDiagnosticDialog()
    % Create modern diagnostic data sharing dialog
    
    % Create figure with better sizing and positioning
    fig = uifigure('Name', 'SCGEATOOL Diagnostic Data', ...
                   'Position', [100 100 450 280], ...
                   'Resize', 'off');
    
    % Create main grid layout with proper spacing
    g = uigridlayout(fig, [5 1]);
    g.RowHeight = {60, '1x', 40, 40, 50};
    g.Padding = [30 25 30 25];
    g.RowSpacing = 15;
    
    % Header panel with icon and title
    headerPanel = uipanel(g, 'BorderType', 'none');
    headerPanel.Layout.Row = 1;
    headerPanel.BackgroundColor = fig.Color;
    
    headerGrid = uigridlayout(headerPanel, [2 1]);
    headerGrid.RowHeight = {'1x', '1x'};
    headerGrid.Padding = [0 0 0 0];
    headerGrid.RowSpacing = 5;
    
    % Title
    titleLabel = uilabel(headerGrid, ...
        'Text', 'Help Improve SCGEATOOL', ...
        'FontSize', 16, ...
        'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center');
    titleLabel.Layout.Row = 1;
    
    % Subtitle
    subtitleLabel = uilabel(headerGrid, ...
        'Text', 'Your privacy matters to us', ...
        'FontSize', 11, ...
        'FontColor', [0.4 0.4 0.4], ...
        'HorizontalAlignment', 'center');
    subtitleLabel.Layout.Row = 2;
    
    % Description panel with wrapped text
    descPanel = uipanel(g, 'BorderType', 'none');
    descPanel.Layout.Row = 2;
    descPanel.BackgroundColor = [0.97 0.97 0.97];
    
    descGrid = uigridlayout(descPanel, [2 2]);
    descGrid.RowHeight = {'1x', 45};
    descGrid.ColumnWidth = {'1x', 60};
    descGrid.Padding = [15 15 15 15];
    descGrid.RowSpacing = 12;
    
    % Main description
    descLabel = uilabel(descGrid, ...
        'Text', ['Share optional diagnostic data to help us improve SCGEATOOL. ' ...
                 'This includes anonymous usage statistics and error reports. ' ...
                 'No personal information is collected.'], ...
        'WordWrap', 'on', ...
        'FontSize', 12);
    descLabel.Layout.Row = 1;
    descLabel.Layout.Column = [1 2];
    
    % Switch label and control in a more prominent layout
    switchLabel = uilabel(descGrid, ...
        'Text', 'Share diagnostic data:', ...
        'FontSize', 12, ...
        'FontWeight', 'bold', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'center');
    switchLabel.Layout.Row = 2;
    switchLabel.Layout.Column = 1;
    
    % uiswitch - the main control element
    diagSwitch = uiswitch(descGrid, 'slider', ...
        'Value', 'On', ...
        'Items', {'Off', 'On'}, ...
        'ValueChangedFcn', @(src, event) switchChanged(src, event));
    diagSwitch.Layout.Row = 2;
    diagSwitch.Layout.Column = 2;
    
    % Checkbox panel
    cbxPanel = uipanel(g, 'BorderType', 'none');
    cbxPanel.Layout.Row = 3;
    cbxPanel.BackgroundColor = fig.Color;
    
    cbxGrid = uigridlayout(cbxPanel, [1 1]);
    cbxGrid.Padding = [0 0 0 0];
    
    dontShowCheckbox = uicheckbox(cbxGrid, ...
        'Text', 'Don''t show this message again', ...
        'Value', false, ...
        'FontSize', 11);
    dontShowCheckbox.Layout.Row = 1;
    dontShowCheckbox.Layout.Column = 1;
    
    % Info label
    infoLabel = uilabel(g, ...
        'Text', 'You can change this preference later in Settings', ...
        'FontSize', 10, ...
        'FontColor', [0.5 0.5 0.5], ...
        'FontAngle', 'italic', ...
        'HorizontalAlignment', 'center');
    infoLabel.Layout.Row = 4;
    
    % Button panel
    btnPanel = uipanel(g, 'BorderType', 'none');
    btnPanel.Layout.Row = 5;
    btnPanel.BackgroundColor = fig.Color;
    
    btnGrid = uigridlayout(btnPanel, [1 3]);
    btnGrid.ColumnWidth = {'1x', 120, 120};
    btnGrid.ColumnSpacing = 10;
    btnGrid.Padding = [0 5 0 0];
    
    % Spacer
    uilabel(btnGrid, 'Text', '');
    
    % Decline button (optional)
    declineBtn = uibutton(btnGrid, ...
        'Text', 'Decline', ...
        'ButtonPushedFcn', @(btn, event) declineCallback(fig, diagSwitch, dontShowCheckbox));
    declineBtn.Layout.Column = 2;
    
    % Continue button
    continueBtn = uibutton(btnGrid, ...
        'Text', 'Continue', ...
        'ButtonPushedFcn', @(btn, event) continueCallback(fig, diagSwitch, dontShowCheckbox), ...
        'BackgroundColor', [0.0 0.4470 0.7410], ...
        'FontColor', [1 1 1], ...
        'FontWeight', 'bold');
    continueBtn.Layout.Column = 3;
    
    % Store UI components in figure for access
    fig.UserData.diagSwitch = diagSwitch;
    fig.UserData.dontShowCheckbox = dontShowCheckbox;
    
    % Center figure on screen
    movegui(fig, 'center');
end

% Callback functions
function switchChanged(src, event)
    if strcmp(src.Value, 'On')
        fprintf('Diagnostic data sharing: ENABLED\n');
    else
        fprintf('Diagnostic data sharing: DISABLED\n');
    end
end

function continueCallback(fig, diagSwitch, dontShowCheckbox)
    % Get user preferences
    shareData = strcmp(diagSwitch.Value, 'On');
    dontShow = dontShowCheckbox.Value;
    
    % Display selection
    fprintf('\n=== User Preferences ===\n');
    fprintf('Share diagnostic data: %s\n', string(shareData));
    fprintf('Don''t show again: %s\n', string(dontShow));
    fprintf('========================\n\n');
    
    % Here you would save preferences to a config file or preferences
    % Example: setpref('SCGEATOOL', 'ShareDiagnostics', shareData);
    %          setpref('SCGEATOOL', 'HideDiagnosticDialog', dontShow);
    
    % Close dialog
    delete(fig);
end

function declineCallback(fig, diagSwitch, dontShowCheckbox)
    % Automatically turn off data sharing and continue
    diagSwitch.Value = 'Off';
    dontShow = dontShowCheckbox.Value;
    
    fprintf('\n=== User Preferences ===\n');
    fprintf('Share diagnostic data: false\n');
    fprintf('Don''t show again: %s\n', string(dontShow));
    fprintf('========================\n\n');
    
    % Close dialog
    delete(fig);
end

% Run the function
createDiagnosticDialog();
%}
