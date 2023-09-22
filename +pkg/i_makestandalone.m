%pw1=fileparts(mfilename('fullpath'));

pw1 = cdgea;
wrkpth = fullfile(pw1, 'tensor_toolbox');
addpath(wrkpth);
wrkpth = fullfile(pw1, '+run', 'external', 'mt_scGeneFit');
addpath(wrkpth);
wrkpth = fullfile(pw1, '+run', 'external', 'mt_MAGIC');
addpath(wrkpth);
wrkpth = fullfile(pw1, '+run', 'external', 'mt_DESeq2');
addpath(wrkpth);
wrkpth = fullfile(pw1, '+run', 'external', 'mt_cbrewer');
addpath(wrkpth);
wrkpth = fullfile(pw1, '+run', 'external', 'mt_PHATE');
addpath(wrkpth);
wrkpth = fullfile(pw1, '+run', 'external', 'mt_UMAP');
addpath(wrkpth);
wrkpth = fullfile(pw1, '+run', 'external', 'mt_SIMLR');
addpath(wrkpth);
wrkpth = fullfile(pw1, '+run', 'external', 'mt_SIMLR', 'src');
addpath(wrkpth);
wrkpth = fullfile(pw1, '+run', 'external', 'mt_CSN_transform');
addpath(wrkpth);
wrkpth = fullfile(pw1, '+run', 'external', 'mt_Specter');
addpath(wrkpth);
wrkpth = fullfile(pw1, '+run', 'external', 'mt_Specter', 'dimred');
addpath(wrkpth);
wrkpth = fullfile(pw1, '+run', 'external', 'mt_Specter', 'LSC');
addpath(wrkpth);
wrkpth = fullfile(pw1, '+run', 'external', 'mt_Specter', 'utils');
addpath(wrkpth);
wrkpth = fullfile(pw1, '+run', 'external', 'mt_SoptSC');
addpath(wrkpth);
wrkpth = fullfile(pw1, '+run', 'external', 'mt_SoptSC', 'NNDSVD');
addpath(wrkpth);
wrkpth = fullfile(pw1, '+run', 'external', 'mt_SoptSC', 'symnmf2');
addpath(wrkpth);
wrkpth = fullfile(pw1, '+run', 'external', 'mt_SinNLRR');
addpath(wrkpth);
wrkpth = fullfile(pw1, '+run', 'external', 'mt_SNNDPC');
addpath(wrkpth);
wrkpth = fullfile(pw1, '+run', 'external', 'mt_SC3', 'ClusterPack');
addpath(wrkpth);
wrkpth = fullfile(pw1, '+run', 'thirdparty', 'locfit', 'm');
addpath(wrkpth);
wrkpth = fullfile(pw1, '+run', 'thirdparty', 'locfit', 'mex');
addpath(wrkpth);
wrkpth = fullfile(pw1, '+run', 'thirdparty', 'alona_panglaodb');
addpath(wrkpth);
wrkpth = fullfile(pw1, '+run', 'thirdparty', 'alona_subtypes');
addpath(wrkpth);


%savepath;

%%
outdir = fullfile(pw1, '..', 'SCGEATOOL_StandaloneApplication');
%
%outdir=sprintf('%s\\Desktop\\scgeatoolstandaloneApplication',a);
if ~exist(outdir, "dir"), mkdir(outdir); end

%%
c = dir('resources/*.gif');
d = strings(length(c), 1);
for k = 1:length(c)
    d(k) = string(fullfile(c(k).folder, c(k).name));
end
%    d=[d;fullfile(pw1, 'resources', 'STRING', 'stringdb_human.mat')];
%    d=[d;fullfile(pw1, 'resources', 'STRING', 'stringdb_mouse.mat')];
d = [d; fullfile(pw1, '+run', 'external', 'mt_UMAP', 'umap.jar')];
d = [d; fullfile(pw1, 'resources', 'tfome_tfgenes.mat')];
d = [d; fullfile(pw1, 'resources', 'cellcyclegenes.xlsx')];
d = [d; fullfile(pw1, 'resources', 'cellscores.xlsx')];
d = [d; fullfile(pw1, 'resources', 'cellscores.txt')];
d = [d; fullfile(pw1, 'resources', 'Ligand_Receptor.mat')];
d = [d; fullfile(pw1, 'resources', 'Ligand_Receptor2.mat')];
d = [d; fullfile(pw1, 'resources', 'myTemplate.pptx')];
d = [d; fullfile(pw1, 'resources', 'DoRothEA_TF_Target_DB', 'dorothea_hs.mat')];
d = [d; fullfile(pw1, 'resources', 'DoRothEA_TF_Target_DB', 'dorothea_mm.mat')];
d = [d; fullfile(pw1, 'example_data', 'testXgs.mat')];
d = [d; fullfile(pw1, 'example_data', 'workshop_example.mat')];
% d=[d;fullfile(matlabroot,'toolbox','rptgen','rptgen','makePPTCompilable.p')];
d2 = string(pkg.dirPlus(fullfile(pw1, '+run', 'external')));
d2 = d2(~contains(d2, "stringdb\stringdb_"));
d = [d; d2];

%%
needcorrect = false;
try
    if ~isdeployed
        compiler.build.standaloneWindowsApplication('scgeatool.m', ...
            'ExecutableName', 'scgeatool', 'Verbose', 'On', ...
            'OutputDir', outdir, 'AdditionalFiles', d, ...
            'SupportPackages', 'autodetect', ...
            'ExecutableVersion', '23.4.4');
    end
catch ME
    disp(ME.message);
    needcorrect = true;
end

%%
try
    a = getenv('USERPROFILE');
    b = getenv('username');

    winopen(sprintf('%s\\AppData\\Local\\Temp\\%s\\mcrCache9.14\\', a, b));
    winopen(outdir);
catch
end

% C:\Users\jcai\AppData\Local\Temp\jcai\mcrCache9.11\scgeat1\scgeatool\+run\external\R_SeuratSaveRds

%%
cd(outdir);
if needcorrect
    a = readmatrix('requiredMCRProducts.txt');
    %writematrix(a(2:end),'requiredMCRProducts.txt','Delimiter','\t');
end
cd ..
zippedfiles = zip('SCGEATOOL_StandaloneApplication.zip', 'SCGEATOOL_StandaloneApplication');
movefile('SCGEATOOL_StandaloneApplication.zip', 'scgeatool.github.io\SCGEATOOL_StandaloneApplication.zip')
