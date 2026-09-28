function [ok, wrkpth, x] = pycommon(prgwkdir)
arguments
prgwkdir{mustBeTextScalar}
end

ok = false;

oldpth = pwd();
cleanupObj = onCleanup(@() cd(oldpth));
pw1 = fileparts(mfilename('fullpath'));
wrkpth = fullfile(pw1, '..',  'external', prgwkdir);
cd(wrkpth);

x = pyenv;
if strlength(x.Executable) == 0, return; end

fw = gui.myWaitbar([], [], [], 'Checking Python environment...');

try
    pkg.i_add_conda_python_path;
catch
    % best-effort: fall back to default pyenv if conda path not found
end

cmdlinestr = sprintf('"%s" "%s%srequire.py"', ...
    x.Executable, wrkpth, filesep);
disp(cmdlinestr)
[status, cmdout] = system(cmdlinestr, '-echo');
if status ~= 0
    if pkg.i_isvalid(fw)
        gui.myWaitbar([], fw, [], 'Checking Python...error.');
    end
    disp(cmdout);
    error('%s has not been installed properly.', ...
        upper(prgwkdir));
end
if pkg.i_isvalid(fw)
    gui.myWaitbar([], fw, [], 'Checking Python environment is complete.');
end
ok = true;
