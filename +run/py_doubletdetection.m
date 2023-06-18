function [isDoublet,doubletscore]=py_doubletdetection(X)
arguments
    X {mustBeNonsparse,mustBeFinite,mustBeNonnegative}
end

isDoublet=[];
doubletscore=[];
    
oldpth=pwd();
[pyok,wrkpth]=run.pycommon('py_doubletdetection');
if ~pyok, return; end


%pw1=fileparts(mfilename('fullpath'));
%wrkpth=fullfile(pw1,'external','py_doubletdetection');
%cd(wrkpth);
%if exist('./input.txt','file'), delete('./input.txt'); end
%if exist('./output1.txt','file'), delete('./output1.txt'); end
%if exist('./output2.txt','file'), delete('./output2.txt'); end

    tmpfilelist={'input.mat','output.mat'};    
    if ~isdebug, pkg.i_deletefiles(tmpfilelist); end


    if issparse(X), X=full(X); end

    save('input.mat','-v7.3','X');
    disp('Input file written.');
    

    cmdlinestr=sprintf('"%s" "%s%sscript.py"', ...
        x.Executable,wrkpth,filesep);
    disp(cmdlinestr)
    [status]=system(cmdlinestr,'-echo');



if status==0 && exist('output.mat','file')
    T=readtable('output1.txt',"ReadVariableNames",false);
    isDoublet=string(T.Var1)=="True"; 
    doubletscore=readmatrix('output2.txt');
end

cd(oldpth);
end

