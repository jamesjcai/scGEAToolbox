function [s]=run_glmpca(X)
%generalized principal component analysis (GLM-PCA)
%
if isempty(FindRpath)
   error('Rscript.exe is not found.');
end
oldpth=pwd;
pw1=fileparts(which(mfilename));
pth=fullfile(pw1,'thirdparty/R_glmpca');
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
%if ~exist('input.csv','file')
    csvwrite('input.csv',transpose(X));
%end
RunRcode('script.R');
if exist('output.csv','file')
    s=csvread('output.csv',1,1);
else
    s=[];
end
if exist('input.csv','file')
    delete('input.csv');
end
if exist('output.csv','file')
    delete('output.csv');
end
cd(oldpth);
