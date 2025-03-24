function [T] = py_memento(wkdir, isdebug)
    T = [];
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
        error(cmdout);
    end
    
    
    %prgfoldername = 'py_writeh5ad';
    %[pyok, wrkpth, x] = run.pycommon(prgfoldername);
    %if ~pyok, return; end
    
    tmpfilelist = {'X.mat', 'g.csv', 'c.csv'};
    if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
    
    
    codefullpath = fullfile(codepth,'script.py');
    pkg.i_addwd2script(codefullpath, wkdir, 'python');
    cmdlinestr = sprintf('"%s" "%s"', x.Executable, codefullpath);
    disp(cmdlinestr)
    [status] = system(cmdlinestr, '-echo');
    
    T = readtable(fullfile(wkdir, 'output.csv'),"FileType","text");
    cd(oldpth);
end