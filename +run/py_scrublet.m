function [isDoublet,doubletscore]=py_scrublet(X)

oldpth=pwd();
isDoublet=[];
doubletscore=[];

    [pyok,wrkpth,x]=run.pycommon('py_scrublet');
    if ~pyok, return; end


%pw1=fileparts(mfilename('fullpath'));
%wrkpth=fullfile(pw1,'external','py_scrublet');
%cd(wrkpth);

    tmpfilelist={'input.mat','output.mat'};
    if ~isdebug, pkg.i_deletefiles(tmpfilelist); end

%if exist('./input.txt','file'), delete('./input.txt'); end
%if exist('./output1.txt','file'), delete('./output1.txt'); end
%if exist('./output2.txt','file'), delete('./output2.txt'); end
if issparse(X), X=full(X); end

%writematrix(X,'input.txt');

%x=pyenv;
%pkg.i_add_conda_python_path;

    save('input.mat','-v7.3','X');
    disp('Input file written.');
    

    cmdlinestr=sprintf('"%s" "%s%sscript.py"', ...
        x.Executable,wrkpth,filesep);
    disp(cmdlinestr)
    [status]=system(cmdlinestr,'-echo');


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

