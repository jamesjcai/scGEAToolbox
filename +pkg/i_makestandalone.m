%pw1=fileparts(mfilename('fullpath'));

pw1=cdgea;
wrkpth=fullfile(pw1,'tensor_toolbox'); addpath(wrkpth);
wrkpth=fullfile(pw1,'+run','thirdparty','cbrewer'); addpath(wrkpth);
wrkpth=fullfile(pw1,'+run','thirdparty','DESeq2'); addpath(wrkpth);
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
wrkpth=fullfile(pw1,'+run','thirdparty','ClusterPack'); addpath(wrkpth);
wrkpth=fullfile(pw1,'+run','thirdparty','MAGIC'); addpath(wrkpth);

%savepath;
%%
outdir = fullfile(pw1, '..','SCGEATOOL_StandaloneApplication');
%
%outdir=sprintf('%s\\Desktop\\scgeatoolstandaloneApplication',a);
if ~exist(outdir,"dir"), mkdir(outdir); end

%%
 c=dir('resources/*.gif');
 d=strings(length(c),1);
 for k=1:length(c)
     d(k)=string(fullfile(c(k).folder,c(k).name));
 end
%    d=[d;fullfile(pw1, 'resources', 'STRING', 'stringdb_human.mat')];
%    d=[d;fullfile(pw1, 'resources', 'STRING', 'stringdb_mouse.mat')];
    d=[d;fullfile(pw1, 'resources', 'tfome_tfgenes.mat')];
    d=[d;fullfile(pw1, 'resources', 'regev_lab_cell_cycle_genes.txt')];
    d=[d;fullfile(pw1, 'resources', 'cellscores.txt')];
    d=[d;fullfile(pw1, 'resources', 'Ligand_Receptor.mat')];
    d=[d;fullfile(pw1, 'resources', 'Ligand_Receptor2.mat')];
    d=[d;fullfile(pw1, 'example_data', 'testSce.mat')];
    d=[d;fullfile(pw1, 'example_data', 'testXgs.mat')];
    d2=string(pkg.dirPlus(fullfile(pw1,'+run','external')));
    d2=d2(~contains(d2,"stringdb\stringdb_"));
    d=[d;d2];

%%
compiler.build.standaloneWindowsApplication('scgeatool.m',...
    'ExecutableName','scgeatool','Verbose','On',...
    'OutputDir',outdir,'AdditionalFiles',d);
%%
try
a=getenv('USERPROFILE');
b=getenv('username');
    
winopen(sprintf('%s\\AppData\\Local\\Temp\\%s\\mcrCache9.11\\',a,b));
winopen(outdir);
catch
end
%%
cd(outdir);
cd ..
zippedfiles = zip('SCGEATOOL_StandaloneApplication.zip','SCGEATOOL_StandaloneApplication');
movefile('SCGEATOOL_StandaloneApplication.zip','scgeatool.github.io\SCGEATOOL_StandaloneApplication.zip')
