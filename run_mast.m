function [T]=run_mast(X,Y,genelist)
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

[~,cmdout]=RunRcode('require.R');
if strfind(cmdout,'there is no package')>0
    cd(oldpth);
    error(cmdout);
end

if exist('output.csv','file')
    delete('output.csv');
end
writematrix(X,'input1.txt');
writematrix(Y,'input2.txt');
RunRcode('script.R');
if exist('output.csv','file')
    T=readtable('output.csv');
    T.Var1=genelist(T.Var1);
    T.Properties.VariableNames{'Var1'} = 'gene';
else
    T=[];
end
if exist('input1.txt','file'), delete('input1.txt'); end
if exist('input2.txt','file'), delete('input2.txt'); end
if exist('output.csv','file'), delete('output.csv'); end
cd(oldpth);
