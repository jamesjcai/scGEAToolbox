function [s]=glmpca(X)
%generalized principal component analysis (GLM-PCA)
%
if isempty(FindRpath)
   error('Rscript.exe is not found.');
end
oldpth=pwd;
pw1=fileparts(mfilename('fullpath'));
pth=fullfile(pw1,'thirdparty','R_glmpca');
cd(pth);
fprintf('CURRENTWDIR = "%s"\n',pth);


[~,cmdout]=pkg.RunRcode('require.R');
if strfind(cmdout,'there is no package')>0
    cd(oldpth);
    error(cmdout);
end

if exist('output.csv','file')
    delete('output.csv');
end
%if ~exist('input.csv','file')
if issparse(X), X=full(X); end
    writematrix(transpose(X),'input.csv');
%end
Rpath=getpref('scgeatoolbox','rexecutablepath');
pkg.RunRcode('script.R',Rpath);

if exist('output.csv','file')
    s=readmatrix('output.csv');
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
end
