%pw1=fileparts(mfilename('fullpath'));
pw1=cdgea;
wrkpth=fullfile(pw1,'+run','thirdparty','cbrewer');
addpath(wrkpth);
wrkpth=fullfile(pw1,'+run','thirdparty','PHATE');
addpath(wrkpth);
wrkpth=fullfile(pw1,'+run','thirdparty','umapFileExchange','umap'); addpath(wrkpth);
wrkpth=fullfile(pw1,'+run','thirdparty','umapFileExchange','util'); addpath(wrkpth);
wrkpth=fullfile(pw1,'+run','thirdparty','R_SeuratWorkflow'); addpath(wrkpth);
wrkpth=fullfile(pw1,'+run','thirdparty','R_MAST'); addpath(wrkpth);
wrkpth=fullfile(pw1,'+run','thirdparty','alona_panglaodb2021'); addpath(wrkpth);

%%
compiler.build.standaloneWindowsApplication('scgeatool.m',...
'ExecutableName','scgeatool')

% C:\Users\jcai\AppData\Local\Temp\jcai\mcrCache9.11

%rootSettings = matlab.internal.getSettingsRoot;
%addonFolder  = rootSettings.matlab.addons.InstallationFolder.ActiveValue;
%executableFolder = fullfile(addonFolder,'Apps',name_of_your_app);
