function [isok] = i_commoncheck_Python(extprogname)
if nargin < 1, extprogname = 'py_GSEApy_enr'; end

    isok = false;
    pw1 = fileparts(mfilename('fullpath'));
    codepth = fullfile(pw1, '..', '+run', 'external', extprogname);
    
    x = pyenv;
    try
        pkg.i_add_conda_python_path;
    catch
    
    end
    
    if ~isempty(extprogname)
        codefullpath = fullfile(codepth,'require.py');
        cmdlinestr = sprintf('"%s" "%s"', x.Executable, codefullpath);
        disp(cmdlinestr)
        [status, cmdout] = system(cmdlinestr, '-echo');
        if status == 0
           isok = true;
        else
            disp(cmdout);
        end
    end
