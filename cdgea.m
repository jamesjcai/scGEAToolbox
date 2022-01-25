function [pw1,pw0]=cdgea(isconfirmed)
%CDGEA - change working directory to the scGEAToolbox folder

if nargin < 1
    isconfirmed=true;
end
pw0=pwd;
pw1=fileparts(mfilename('fullpath'));



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

if ~exist(fullfile(pw1,'example_data/'),'dir') || ~exist(fullfile(pw1,'resources/'),'dir')
    f = waitbar(0.33,'Initializing scgeatool on the first run...');
    if ~exist(fullfile(pw1,'example_data/'),'dir')
        try
            unzip('https://github.com/jamesjcai/jamesjcai.github.io/raw/master/data/example_data.zip');
        catch ME
            close(f);
            warndlg(ME.message);
            warning(ME.message);
        end
    end
    f = waitbar(0.66,'Initializing scgeatool on the first run...');
    if ~exist(fullfile(pw1,'resources/'),'dir')
        try
            unzip('https://github.com/jamesjcai/jamesjcai.github.io/raw/master/data/resources.zip');
        catch ME
            close(f);
            warndlg(ME.message);
            warning(ME.message);
        end
    end
    waitbar(1,f,'Finishing');
    close(f);
end
% if ~exist(fullfile(pw1,'+run','thirdparty/'),'dir')
%     try
%     unzip('https://github.com/jamesjcai/jamesjcai.github.io/raw/master/data/thirdparty.zip');
%     movefile(fullfile(pw1,'thirdparty/'),fullfile(pw1,'+run','thirdparty/'));
%     catch ME
%         warndlg(ME.message);
%         warning(ME.message);
%     end
% end
end
