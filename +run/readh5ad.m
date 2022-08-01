function [X,g,b]=readh5ad(fname)


if ~exist(fname,'file')
    error('FNAME is not a file.');
end

oldpth=pwd();
pw1=fileparts(mfilename('fullpath'));
wrkpth=fullfile(pw1,'external','py_readh5ad');
cd(wrkpth);
tmpfilelist={'output.mat'};
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end

x=pyenv;
pkg.i_add_conda_python_path;
cmdlinestr=sprintf('"%s" "%s%sscript_csv.py" "%s"',x.Executable, ...
        wrkpth,filesep,fname);
disp(cmdlinestr)
[status]=system(cmdlinestr);

if status==0 && exist('output.txt','file')

end

if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
cd(oldpth);
end
