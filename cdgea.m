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

% if ~exist(fullfile(pw1,'example_data/'),'dir') || ~exist(fullfile(pw1,'resources/'),'dir') ||...
%         ~exist(fullfile(pw1,'+run','external/'),'dir')
%     f = waitbar(0.33,'Initializing SCGEATOOL on the first run...');
%     if ~exist(fullfile(pw1,'example_data/'),'dir')
%         try
%             unzip('https://github.com/jamesjcai/jamesjcai.github.io/raw/master/data/example_data.zip');
%         catch ME
%             warndlg(ME.message);
%             warning(ME.message);
%         end
%     end
%     waitbar(0.66,f,'Initializing SCGEATOOL on the first run...');
%     if ~exist(fullfile(pw1,'resources/'),'dir')
%         try
%             unzip('https://github.com/jamesjcai/jamesjcai.github.io/raw/master/data/resources.zip');
%         catch ME
%             warndlg(ME.message);
%             warning(ME.message);
%         end
%     end
%     waitbar(0.85,f,'Initializing SCGEATOOL on the first run...');
%     if ~exist(fullfile(pw1,'+run','external/'),'dir')
%         try
%             a=unzip('https://github.com/jamesjcai/jamesjcai.github.io/raw/master/data/external.zip');
%         catch ME
%             warndlg(ME.message);
%             warning(ME.message);
%         end            
%             if ~iscell(a)
%                 errordlg('SCGEATOOL could not download the supporting file.\nPlease check network connection and try again.');
%                 return;
%             end
%         try
%             movefile(fullfile(pw1,'external/'),fullfile(pw1,'+run','external/'));
%         catch ME
%             warndlg(ME.message);
%             warning(ME.message);
%         end
%     end
%     pause(0.5);
%     waitbar(1,f,'Finishing');
%     pause(1);
%     if isvalid(f), close(f); end
% end

end
