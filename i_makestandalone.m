pw1=fileparts(mfilename('fullpath'));
wrkpth=fullfile(pw1,'+run','thirdparty','cbrewer');
addpath(wrkpth);
wrkpth=fullfile(pw1,'+run','thirdparty','PHATE');
addpath(wrkpth);
wrkpth=fullfile(pw1,'+run','thirdparty','umapFileExchange');
addpath(wrkpth);
wrkpth=fullfile(pw1,'+run','thirdparty','R_SeuratWorkflow');
addpath(wrkpth);
wrkpth=fullfile(pw1,'+run','thirdparty','R_MAST');
addpath(wrkpth);
wrkpth=fullfile(pw1,'+run','thirdparty','alona_panglaodb2021');
addpath(wrkpth);


compiler.build.standaloneWindowsApplication('scgeatool.m',...
'ExecutableName','scgeatool')

% C:\Users\jcai\AppData\Local\Temp\jcai\mcrCache9.11