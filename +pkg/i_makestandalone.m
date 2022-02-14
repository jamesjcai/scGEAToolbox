%pw1=fileparts(mfilename('fullpath'));

pw1=cdgea;
wrkpth=fullfile(pw1,'+run','thirdparty','tensor_toolbox'); addpath(wrkpth);
wrkpth=fullfile(pw1,'+run','thirdparty','cbrewer'); addpath(wrkpth);
wrkpth=fullfile(pw1,'+run','thirdparty','PHATE'); addpath(wrkpth);
wrkpth=fullfile(pw1,'+run','thirdparty','umapFileExchange'); addpath(wrkpth);
wrkpth=fullfile(pw1,'+run','thirdparty','alona_panglaodb'); addpath(wrkpth);
wrkpth=fullfile(pw1,'+run','thirdparty','alona_subtypes'); addpath(wrkpth);
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

%%
 c=dir('resources/*.gif');
 d=strings(length(c),1);
 for k=1:length(c)
     d(k)=string(fullfile(c(k).folder,c(k).name));
 end
    d=[d;fullfile(pw1, 'resources', 'STRING', 'stringdb_human.mat')];
    d=[d;fullfile(pw1, 'resources', 'STRING', 'stringdb_mouse.mat')];
    d=[d;fullfile(pw1, 'resources', 'tfome_tfgenes.mat')];
    d=[d;fullfile(pw1, 'resources', 'regev_lab_cell_cycle_genes.txt')];
    d=[d;fullfile(pw1, 'resources', 'cellscores.txt')];
    d=[d;fullfile(pw1, 'resources', 'Ligand_Receptor.mat')];
    d=[d;fullfile(pw1, 'resources', 'Ligand_Receptor2.mat')];
    d=[d;fullfile(pw1, 'example_data', 'testSce.mat')];
    d=[d;fullfile(pw1, 'example_data', 'testXgs.mat')];

%%
compiler.build.standaloneWindowsApplication('scgeatool.m',...
    'ExecutableName','scgeatool','Verbose','On',...
    'OutputDir',outdir,'AdditionalFiles',d);
%%
winopen(sprintf('%s\\AppData\\Local\\Temp\\%s\\mcrCache9.11\\',a,b));
winopen(outdir);




