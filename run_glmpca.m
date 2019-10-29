function [s]=run_glmpca(X)
if isempty(FindRpath)
   error('Rscript.exe is not found.');
end
oldpth=pwd;
pw1=fileparts(which(mfilename));
pth=fullfile(pw1,'thirdparty/R_glmpca');
cd(pth);

if exist('output.csv','file')
    delete('output.csv');
end
%if ~exist('input.csv','file')
    csvwrite('input.csv',X);
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
