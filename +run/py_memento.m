function [succeeded] = py_memento(wkdir, isdebug)

succeeded = false;
extprogname = 'py_memento';
if nargin<1 || isempty(wkdir)
    preftagname = 'externalwrkpath';
    [wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname);
    if isempty(wkdir), return; end
end
if nargin < 2, isdebug = true; end

oldpth = pwd();
pw1 = fileparts(mfilename('fullpath'));
codepth = fullfile(pw1, 'external', extprogname);
if isempty(wkdir) || ~isfolder(wkdir)
    cd(codepth);
else
    disp('Using working directory provided.');
    cd(wkdir);
end
fw = gui.gui_waitbar([], [], 'Checking Python environment...');

x = pyenv;
try
    pkg.i_add_conda_python_path;
catch

end

codefullpath = fullfile(codepth,'require.py');
cmdlinestr = sprintf('"%s" "%s"', x.Executable, codefullpath);
disp(cmdlinestr)
[status, cmdout] = system(cmdlinestr, '-echo');
if status ~= 0
    cd(oldpth);
    if isvalid(fw)
        gui.gui_waitbar(fw, true);
    end
    error(cmdout);
end


%prgfoldername = 'py_writeh5ad';
%[pyok, wrkpth, x] = run.pycommon(prgfoldername);
%if ~pyok, return; end

tmpfilelist = {'X.mat', 'g.csv', 'c.csv'};
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end


if isvalid(fw)
    gui.gui_waitbar(fw, [], [], 'Checking Python environment is complete');
    pause(0.5);
    gui.gui_waitbar(fw, [], [], 'Running py\_memento...');
end

codefullpath = fullfile(codepth,'script.py');
pkg.i_addwd2script(codefullpath, wkdir, 'python');
cmdlinestr = sprintf('"%s" "%s"', x.Executable, codefullpath);
disp(cmdlinestr)
[status] = system(cmdlinestr, '-echo');

if status == 0 && isvalid(fw)
    gui.gui_waitbar(fw, false, 'File is written.');
    succeeded = true;
else
    gui.gui_waitbar(fw, true, 'File is failed to save.');
end

cd(oldpth);
end