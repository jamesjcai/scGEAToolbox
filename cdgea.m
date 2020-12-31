function [pw1,pw0]=cdgea(isconfirmed)
%CDGEA Changes to GEAToolbox directory

if nargin < 1
    isconfirmed=true;
end
pw0=pwd;
pw1=fileparts(which(mfilename));
if ~strcmp(pw0,pw1) && ~isconfirmed
    [selectedButton]=uigetpref('scGEApp',... % Group
           'cdgea_ask',...                               % Preference
           'Changing Working Directory',...              % Window title
           {'Do you want to change current working directory to scGEApp directory?'},...
           {'always','never';'Yes','No'},...       % Values and button strings
           'ExtraOptions','Cancel',...             % Additional button
           'DefaultButton','Yes');
    switch selectedButton
        case {'always','Yes'}
            cd(pw1);
        case {'never','No','Cancel'}
            % do nothing
    end
else
    cd(pw1);
end