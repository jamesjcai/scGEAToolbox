function [isDoublet]=doubletdetection(X)

oldpth=pwd();
pw1=fileparts(mfilename('fullpath'));
wrkpth=fullfile(pw1,'thirdparty','doubletdetection');
cd(wrkpth);

if exist('input.txt','file'), delete('input.txt'); end
if exist('output.txt','file'), delete('output.txt'); end

writematrix(X,'input.txt');


x=pyenv;
pkg.i_add_conda_python_path;
cmdlinestr=sprintf('"%s" "%s%sscript.py"',x.Executable,wrkpth,filesep);
disp(cmdlinestr)
[status]=system(cmdlinestr);

if status==0 && exist('output.txt','file')
    T=readtable('output.txt',"ReadVariableNames",false);
    isDoublet=string(T.Var1)=="True";    
else
    isDoublet=[];
end

if exist('input.txt','file'), delete('input.txt'); end
if exist('output.txt','file'), delete('output.txt'); end
cd(oldpth);
end

