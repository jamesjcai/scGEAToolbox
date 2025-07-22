%pw1=fileparts(mfilename('fullpath'));

pw1 = cdgea;
wrkpth = fullfile(pw1, 'external', 'tensor_toolbox');
addpath(wrkpth);
wrkpth = fullfile(pw1, 'external', 'ml_scGeneFit');
addpath(wrkpth);
wrkpth = fullfile(pw1, 'external', 'ml_MAGIC');
addpath(wrkpth);
wrkpth = fullfile(pw1, 'external', 'ml_cbrewer');
addpath(wrkpth);
wrkpth = fullfile(pw1, 'external', 'ml_PHATE');
addpath(wrkpth);
wrkpth = fullfile(pw1, 'external', 'ml_UMAP44');
addpath(wrkpth);
wrkpth = fullfile(pw1, 'external', 'ml_SIMLR');
addpath(wrkpth);
wrkpth = fullfile(pw1, 'external', 'ml_SIMLR', 'src');
addpath(wrkpth);
wrkpth = fullfile(pw1, 'external', 'ml_CSN_transform');
addpath(wrkpth);
wrkpth = fullfile(pw1, 'external', 'ml_Specter');
addpath(wrkpth);
wrkpth = fullfile(pw1, 'external', 'ml_Specter', 'dimred');
addpath(wrkpth);
wrkpth = fullfile(pw1, 'external', 'ml_Specter', 'LSC');
addpath(wrkpth);
wrkpth = fullfile(pw1, 'external', 'ml_Specter', 'utils');
addpath(wrkpth);
wrkpth = fullfile(pw1, 'external', 'ml_SoptSC');
addpath(wrkpth);
wrkpth = fullfile(pw1, 'external', 'ml_SoptSC', 'NNDSVD');
addpath(wrkpth);
wrkpth = fullfile(pw1, 'external', 'ml_SoptSC', 'symnmf2');
addpath(wrkpth);
wrkpth = fullfile(pw1, 'external', 'ml_SinNLRR');
addpath(wrkpth);
wrkpth = fullfile(pw1, 'external', 'ml_SNNDPC');
addpath(wrkpth);
wrkpth = fullfile(pw1, 'external', 'ml_SC3', 'ClusterPack');
addpath(wrkpth);
wrkpth = fullfile(pw1, 'external', 'locfit');
addpath(wrkpth);
wrkpth = fullfile(pw1, 'external', 'alona_panglaodb');
addpath(wrkpth);
wrkpth = fullfile(pw1, 'external', 'alona_subtypes');
addpath(wrkpth);


%savepath;

%%
outdir = fullfile(pw1, '..', 'SCGEATOOL_StandaloneApplication');
%
%outdir=sprintf('%s\\Desktop\\scgeatoolstandaloneApplication',a);
if ~exist(outdir, "dir"), mkdir(outdir); end

%%
c = dir('assets/Images/*.gif');
d1 = strings(length(c), 1);
for k = 1:length(c)
    d1(k) = string(fullfile(c(k).folder, c(k).name));
end
c = dir('assets/Images/*.png');
d2 = strings(length(c), 1);
for k = 1:length(c)
    d2(k) = string(fullfile(c(k).folder, c(k).name));
end
c = dir('assets/Images/*.jpg');
d3 = strings(length(c), 1);
for k = 1:length(c)
    d3(k) = string(fullfile(c(k).folder, c(k).name));
end
c = dir('assets/Images/*.mat');
d4 = strings(length(c), 1);
for k = 1:length(c)
    d4(k) = string(fullfile(c(k).folder, c(k).name));
end
d = [d1;d2;d3;d4];
%    d=[d;fullfile(pw1, 'assets', 'STRING', 'stringdb_human.mat')];
%    d=[d;fullfile(pw1, 'assets', 'STRING', 'stringdb_mouse.mat')];
% d = [d; fullfile(pw1, 'external', 'ml_UMAP', 'umap.jar')];
d = [d; fullfile(pw1, 'assets', 'Misc', 'refinfo.txt')];
d = [d; fullfile(pw1, 'assets', 'TFome', 'tfome_tfgenes.mat')];
d = [d; fullfile(pw1, 'assets', 'CellScores', 'cellcyclegenes.xlsx')];
d = [d; fullfile(pw1, 'assets', 'CellScores', 'cellscores.xlsx')];
d = [d; fullfile(pw1, 'assets', 'CellScores', 'cellscores.txt')];
d = [d; fullfile(pw1, 'assets', 'ScTypeDB', 'ScTypeDB_full.xlsx')];
d = [d; fullfile(pw1, 'assets', 'PanglaoDB', 'celltypes.xlsx')];
%d = [d; fullfile(pw1, 'assets', 'Ligand_Receptor.mat')];
%d = [d; fullfile(pw1, 'assets', 'Ligand_Receptor2.mat')];
d = [d; fullfile(pw1, 'assets', 'Ligand_Receptor','Ligand_Receptor_more.mat')];
d = [d; fullfile(pw1, 'assets', 'Misc', 'myTemplate.pptx')];
d = [d; fullfile(pw1, 'assets', 'Misc', 'myTemplate2.pptx')];
d = [d; fullfile(pw1, 'assets', 'Misc', 'value_template_pos.txt')];
d = [d; fullfile(pw1, 'assets', 'Misc', 'value_template_std.txt')];
d = [d; fullfile(pw1, 'assets', 'DoRothEA_TF_Target_DB', 'dorothea_hs.mat')];
d = [d; fullfile(pw1, 'assets', 'DoRothEA_TF_Target_DB', 'dorothea_mm.mat')];
d = [d; fullfile(pw1, 'example_data', 'testXgs.mat')];
d = [d; fullfile(pw1, 'example_data', 'new_example_sce.mat')];
d = [d; fullfile(pw1, 'scGEAToolbox.prj')];

% d=[d;fullfile(matlabroot,'toolbox','rptgen','rptgen','makePPTCompilable.p')];
d2 = string(dirPlus(fullfile(pw1, 'external')));
d2 = d2(~contains(d2, "stringdb\stringdb_"));
d = [d; d2];

%%
needcorrect = false;
try
    if ~isdeployed
        % a = datestr(datetime("today"), 'yy.mm.dd');
        a = string(datetime("today",'Format','yy.MM.dd'));
        compiler.build.standaloneWindowsApplication('scgeatool.m', ...
            'ExecutableName', 'scgeatool', 'Verbose', 'On', ...
            'OutputDir', outdir, 'AdditionalFiles', d, ...
            'SupportPackages', 'autodetect', ...
            'ExecutableVersion', a);
    end
catch ME
    disp(ME.message);
    needcorrect = true;
end

%%
try
    a = getenv('USERPROFILE');
    b = getenv('username');
    winopen(sprintf('%s\\AppData\\Local\\Temp\\%s\\mcrCache23.2\\', a, b));
    winopen(outdir);
catch
end

% C:\Users\jcai\AppData\Local\Temp\jcai\mcrCache9.11\scgeat1\scgeatool\+run\external\R_SeuratSaveRds

%%

copyfile('assets/Images/splash.png', fullfile(outdir,"splash.png"));
cd(outdir);
% if needcorrect
%     a = readmatrix('requiredMCRProducts.txt');
%     %writematrix(a(2:end),'requiredMCRProducts.txt','Delimiter','\t');
% end
cd ..
zippedfiles = zip('SCGEATOOL_StandaloneApplication.zip', 'SCGEATOOL_StandaloneApplication');
movefile('SCGEATOOL_StandaloneApplication.zip', 'scgeatool.github.io\SCGEATOOL_StandaloneApplication.zip')
cdgea;
