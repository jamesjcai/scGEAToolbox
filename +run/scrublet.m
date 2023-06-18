function [isDoublet,doubletscore]=scrublet(X)

oldpth=pwd();
pw1=fileparts(mfilename('fullpath'));
wrkpth=fullfile(pw1,'external','py_scrublet');
cd(wrkpth);

if exist('./input.txt','file'), delete('./input.txt'); end
if exist('./output1.txt','file'), delete('./output1.txt'); end
if exist('./output2.txt','file'), delete('./output2.txt'); end
if issparse(X), X=full(X); end
%writematrix(X,'input.txt');
    save('input.mat','-v7.3','X');

x=pyenv;
pkg.i_add_conda_python_path;
cmdlinestr=sprintf('"%s" "%s%sscript.py"',x.Executable,wrkpth,filesep);
disp(cmdlinestr)
[status]=system(cmdlinestr);

if status==0 && exist('output1.txt','file')
    T=readtable('output1.txt',"ReadVariableNames",false);
    isDoublet=string(T.Var1)=="True"; 
    doubletscore=readmatrix('output2.txt');
else
    isDoublet=[];
    doubletscore=[];
end

if exist('./input.txt','file'), delete('./input.txt'); end
if exist('./output1.txt','file'), delete('./output1.txt'); end
if exist('./output2.txt','file'), delete('./output2.txt'); end
cd(oldpth);
end

