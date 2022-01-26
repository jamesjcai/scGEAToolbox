%pw1=fileparts(mfilename('fullpath'));

pw1=cdgea;
wrkpth=fullfile(pw1,'+run','thirdparty','cbrewer'); addpath(wrkpth);
wrkpth=fullfile(pw1,'+run','thirdparty','PHATE'); addpath(wrkpth);
wrkpth=fullfile(pw1,'+run','thirdparty','umapFileExchange'); addpath(wrkpth);
wrkpth=fullfile(pw1,'+run','thirdparty','alona_panglaodb2021'); addpath(wrkpth);
wrkpth=fullfile(pw1,'+run','thirdparty','SIMLR'); addpath(wrkpth);
wrkpth=fullfile(pw1,'+run','thirdparty','SIMLR','src'); addpath(wrkpth);
wrkpth=fullfile(pw1,'+run','thirdparty','CSN_transform'); addpath(wrkpth);
wrkpth=fullfile(pw1,'+run','thirdparty','scGeneFit'); addpath(wrkpth);
wrkpth=fullfile(pw1,'+run','thirdparty','Specter'); addpath(wrkpth);
wrkpth=fullfile(pw1,'+run','thirdparty','Specter','dimred'); addpath(wrkpth);
wrkpth=fullfile(pw1,'+run','thirdparty','Specter','LSC'); addpath(wrkpth);
wrkpth=fullfile(pw1,'+run','thirdparty','Specter','utils'); addpath(wrkpth);
wrkpth=fullfile(pw1,'+run','thirdparty','SoptSC'); addpath(wrkpth);
wrkpth=fullfile(pw1,'+run','thirdparty','SoptSC','NNDSVD'); addpath(wrkpth);
wrkpth=fullfile(pw1,'+run','thirdparty','SoptSC','symnmf2'); addpath(wrkpth);
wrkpth=fullfile(pw1,'+run','thirdparty','SNNDPC'); addpath(wrkpth);
wrkpth=fullfile(pw1,'+run','thirdparty','SinNLRR'); addpath(wrkpth);

%savepath;
%%
compiler.build.standaloneWindowsApplication('scgeatool.m',...
    'ExecutableName','scgeatool','Verbose','On');

% C:\Users\jcai\AppData\Local\Temp\jcai\mcrCache9.11

%rootSettings = matlab.internal.getSettingsRoot;
%addonFolder  = rootSettings.matlab.addons.InstallationFolder.ActiveValue;
%executableFolder = fullfile(addonFolder,'Apps',name_of_your_app);
