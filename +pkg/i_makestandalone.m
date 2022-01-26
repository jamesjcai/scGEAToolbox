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
a=getenv('USERPROFILE');
b=getenv('username');
outdir=sprintf('%s\\Desktop\\scgeatoolstandaloneApplication',a);
if ~exist(outdir,"dir"), makedir(outdir); end

compiler.build.standaloneWindowsApplication('scgeatool.m',...
    'ExecutableName','scgeatool','Verbose','On',...
    'OutputDir',outdir);
%%
winopen(sprintf('%s\\AppData\\Local\\Temp\\%s\\mcrCache9.11\\',a,b));
winopen(outdir);


