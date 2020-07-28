function [A]=run_mast(X,Y)

% if nargin<2, plotit=false; end
if isempty(FindRpath)
   error('Rscript.exe is not found.');
end

oldpth=pwd;
pw1=fileparts(which(mfilename));
pth=fullfile(pw1,'thirdparty/R_MAST');
cd(pth);
fprintf('CURRENTWDIR = "%s"\n',pth);

writematrix(X,'input1.txt');
writematrix(Y,'input2.txt');
% RunRcode('script.R');
if exist('output.txt','file')
    A=readmatrix('output.txt',1,1);
else
    A=[];
end
cd(oldpth);
