function [celltypes, T] = py_panhumanpy(sce, wkdir, ...
    isdebug, prepare_input_only)

% cell_type = run.py_panhumanpy(sce, 'C:\Users\jcai\Downloads', true, true);

celltypes = [];
if nargin < 4, prepare_input_only = false; end
if nargin < 3, isdebug = true; end
if nargin < 2, wkdir = tempdir; end

oldpth = pwd();
pw1 = fileparts(mfilename('fullpath'));
codepth = fullfile(pw1, '..', 'external', 'py_panhumanpy');

if isempty(wkdir) || ~isfolder(wkdir)
    cd(codepth);
else
    disp('Using working directory provided.');
    cd(wkdir);
end
% winopen(wkdir);

% fw = gui.gui_waitbar([], [], 'Checking Python environment...');

x = pyenv;
try
    pkg.i_add_conda_python_path;
catch

end

codepth = pkg.i_normalizepath(codepth);

    if ~prepare_input_only

        codefullpath = fullfile(codepth,'require.py');
        cmdlinestr = sprintf('"%s" "%s"', x.Executable, codefullpath);
        
        disp(cmdlinestr)
        [status, cmdout] = system(cmdlinestr, '-echo');
        if status ~= 0
            cd(oldpth);
            error(cmdout);
        else 
            disp('Code requirement check is done.')
        end
    end

%try
    pkg.i_deletefiles({'input.h5ad', 'output.h5ad','tg.csv'});
    tmpfilelist = {'Xnorm.mat', 'X.mat', 'g.csv', 'c.csv', 'tg.csv', ...
        'input.h5ad', 'output.h5ad'};
    if ~isdebug, pkg.i_deletefiles(tmpfilelist); end

    [Xnorm] = pkg.norm_libsize(sce.X, 10000);
    Xnorm = log1p(full(Xnorm));
    Xnorm = single(Xnorm);

    g = cellstr(sce.g);
    save('X.mat','-v7.3',"Xnorm","g");
    
codefullpath = fullfile(codepth,'script.py');
pkg.i_addwd2script(codefullpath, wkdir, 'python');
cmdlinestr = sprintf('"%s" "%s"', x.Executable, codefullpath);
disp(cmdlinestr)

if ~prepare_input_only
    [status] = system(cmdlinestr, '-echo');
    if status == 0 && exist('output.csv', 'file')
        T = readtable('output.csv','ReadVariableNames', true, ...
            'VariableNamingRule', 'modify');
        celltypes = string(T.final_level_labels);
    end
else
    disp('Input files are prepared. To do the analysis, run script.py in the working folder.')
end

if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
cd(oldpth);
end
