function [c]=run_cellcycle(X,genelist)
%Run cell cycle analysis using R package Seurat

if nargin<2, error("[c]=run_cellcycle(X,genelist)"); end

if isempty(FindRpath)
   error('Rscript.exe is not found.');
end
oldpth=pwd;
pw1=fileparts(which(mfilename));
pth=fullfile(pw1,'thirdparty/R_cellcycle');
cd(pth);
fprintf('CURRENTWDIR = "%s"\n',pth);


[~,cmdout]=RunRcode('require.R');
if strfind(cmdout,'there is no package')>0
    cd(oldpth);
    error(cmdout);
end

if exist('output.csv','file')
    delete('output.csv');
end
sc_writefile('input.txt',X,upper(genelist));
RunRcode('script.R');


if exist('output.csv','file')
    T=readtable('output.csv','ReadVariableNames',true);
    c=string(T.Phase);
else
    c=[];
end
if exist('input.txt','file'), delete('input.txt'); end
if exist('output.csv','file'), delete('output.csv'); end
cd(oldpth);
end