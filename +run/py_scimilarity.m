function [s] = py_scimilarity(sce, wkdir, isdebug)
s = [];
if nargin < 3, isdebug = true; end
if nargin < 2, wkdir = tempdir; end

oldpth = pwd();
pw1 = fileparts(mfilename('fullpath'));
codepth = fullfile(pw1, 'external', 'py_scimilarity');

if isempty(wkdir) || ~isfolder(wkdir)
    cd(codepth);
else
    disp('Using working directory provided.');
    cd(wkdir);
end
winopen(wkdir);

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
else 
    disp('Code requirement check is done.')
end

try
    pkg.i_deletefiles({'input.h5ad', 'output.h5ad'});
    tmpfilelist = {'Xnorm.mat', 'X.mat', 'g.csv', 'c.csv', ...
        'input.h5ad', 'output.h5ad'};
    if ~isdebug, pkg.i_deletefiles(tmpfilelist); end    
    if issparse(sce.X)
        X = single(full(sce.X));         
    else
        X = single(sce.X);
    end
    save('X.mat', '-v7.3', 'X');
    Xnorm = single(sc_norm(full(sce.X)));
    save('Xnorm.mat','-v7.3',"Xnorm");
    g = sce.g;
    writetable(table(g),'g.csv','WriteVariableNames',false);
    sce.c_cell_id = matlab.lang.makeUniqueStrings(sce.c_cell_id);
    T = pkg.makeattributestable(sce);
    writetable(T,'c.csv');
catch ME
    if isvalid(fw)
         gui.gui_waitbar(fw, true);
    end
    errordlg(ME.message,'');
    return;
end
if isvalid(fw)
    gui.gui_waitbar(fw, [], [], 'Checking Python environment is complete');
    pause(0.5);
    gui.gui_waitbar(fw, [], [], sprintf('Running %s...', 'py_scimilarity'));
end

codefullpath = fullfile(codepth,'script.py');
pkg.i_addwd2script(codefullpath, wkdir, 'python');
cmdlinestr = sprintf('"%s" "%s"', x.Executable, codefullpath);
disp(cmdlinestr)
[status] = system(cmdlinestr, '-echo');
% [status2] = movefile('output.h5ad',fname);
if status == 0 && isvalid(fw)
    gui.gui_waitbar(fw, [], 'output.h5ad is written.');
end
if status == 0 && exist('output.h5ad', 'file')
    % adams.obs["predictions_unconstrained"]
    s = h5read('output.h5ad','/obs/predictions_unconstrained');
    sce.c_celltype_tx = s;
end

if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
cd(oldpth);

end