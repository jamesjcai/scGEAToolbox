function [regNet, tfgenes, targetg] = r_SCORPION(X, g, wkdir, isdebug)
% Run SCORPION
if nargin < 3, wkdir = tempdir; end
if nargin < 4, isdebug = true; end

tfgenes = [];
targetg = [];
regNet = [];
oldpth = pwd();
[isok, msg, codepth] = commoncheck_R('R_SCORPION');
if ~isok, error(msg); end
if ~isempty(wkdir) && isfolder(wkdir), cd(wkdir); end

tmpfilelist = {'input.h5', 'output.h5', 'hg38_PPI.RData','hg38_TF.RData'};
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end

%pw1 = fileparts(mfilename('fullpath'));
%fname = fullfile(pw1, '..', 'resources', 'DoRothEA_TF_Target_DB', 'dorothea_hs.mat');
%load(fname,'T');
%writetable(T(:,{'tf','target','mor'}),'tf.csv');
%dbfile1 = fullfile(pw1, 'external', 'web_STRING', 'stringdb', 'stringdb_human.mat');
%load(dbfile1,'G');
%G = G.subgraph(ismember(G.Nodes.Name,g));


websave('hg38_PPI.RData','https://github.com/dosorio/SCORPION/raw/main/Data/hg38_PPI.RData');
websave('hg38_TF.RData','https://github.com/dosorio/SCORPION/raw/main/Data/hg38_TF.RData');

if issparse(X), X = full(X); end
pkg.e_writeh5(X, g, 'input.h5');

Rpath = getpref('scgeatoolbox', 'rexecutablepath',[]);
if isempty(Rpath)
    error('R environment has not been set up.');
end

codefullpath = fullfile(codepth,'script.R');
pkg.RunRcode(codefullpath, Rpath);
if exist('output.h5', 'file')
    regNet = h5read('output.h5', '/regNet');
    %tfgenes =  strrep(h5read('output.h5', '/rownames'),char(0),"");
    %targetg = strrep(h5read('output.h5', '/colnames'),char(0),"");
    tfgenes = deblank(h5read('output.h5', '/rownames'));
    targetg = deblank(h5read('output.h5', '/colnames'));

    T = array2table(regNet);
    T.Row=tfgenes;
    T.Properties.VariableNames=targetg;
    writetable(T,'output.txt','WriteRowNames',true);

    %[a,b]=maxk(regNet(:),5);
    %[x,y]=ind2sub(size(regNet),b);
    %[tfgenes(x) targetg(y)]

end
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
cd(oldpth);
end
