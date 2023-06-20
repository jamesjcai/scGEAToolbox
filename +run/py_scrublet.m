function [isDoublet,doubletscore]=py_scrublet(X)

oldpth=pwd();
isDoublet=[];
doubletscore=[];
prgfoldername='py_scrublet';
[pyok,wrkpth,x]=run.pycommon(prgfoldername);
if ~pyok, return; end
isdebug=false;

tmpfilelist={'input.mat','output.mat'};
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end

if issparse(X), X=full(X); end

save('input.mat','-v7.3','X');
disp('Input file written.');
    
% cmdlinestr=sprintf('"%s" "%s%sscript.py"', ...
%     x.Executable,wrkpth,filesep);
% disp(cmdlinestr)
% [status]=system(cmdlinestr,'-echo');
    [status]=run.pycommon2(x,wrkpth,prgfoldername);

if status==0 && exist('output.mat','file')
    load("output.mat",'isDoublet','doubletscore')
end
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
cd(oldpth);
end
