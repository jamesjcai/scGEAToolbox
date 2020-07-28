function [A]=run_mast(X,Y,genelist)
if nargin<3
    genelist=(1:size(X,1))';
end
if isempty(FindRpath)
   error('Rscript.exe is not found.');
end
oldpth=pwd;
pw1=fileparts(which(mfilename));
pth=fullfile(pw1,'thirdparty/R_MAST');
cd(pth);
fprintf('CURRENTWDIR = "%s"\n',pth);
if exist('output.csv','file')
    delete('output.csv');
end
writematrix(X,'input1.txt');
writematrix(Y,'input2.txt');
% RunRcode('script.R');
if exist('output.csv','file')
    A=readtable('output.csv');
    A.Var1=genelist(A.Var1);
    A.Properties.VariableNames{'Var1'} = 'gene';
else
    A=[];
end
cd(oldpth);
