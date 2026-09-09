function [succeeded] = py_writeh5ad(sce, fname, wkdir, isdebug, verbose, pe)

if nargin < 6, pe = []; end

succeeded = false;
if nargin < 2, fname = tempname + ".h5ad"; end
extprogname = 'py_writeh5ad';
if nargin<3 || isempty(wkdir)
    preftagname = 'externalwrkpath';
    [wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname);
    if isempty(wkdir), return; end
end
if nargin < 4, isdebug = true; end
if nargin < 5, verbose = true; end

oldpth = pwd();
cleanupCwd = onCleanup(@() cd(oldpth));
pw1 = fileparts(mfilename('fullpath'));
codepth = fullfile(pw1, '..',  'external', extprogname);
if isempty(wkdir) || ~isfolder(wkdir)
    cd(codepth);
else
    if verbose
        disp('Using working directory provided.');
    end
    cd(wkdir);
end

if verbose
    fw = gui.gui_waitbar([], [], 'Checking Python environment...');
end

if isempty(pe)
    pe = pyenv;
end

% pyenv( ...
%     "Version","D:\claude_dev\GEOcellar_autogen\.venv\Scripts\python.exe", ...
%     "ExecutionMode","InProcess")

% try
%     pkg.i_add_conda_python_path;
% catch
% 
% end
codepth = pkg.i_normalizepath(codepth);

if verbose
    codefullpath = fullfile(codepth,'require.py');
    cmdlinestr = sprintf('"%s" "%s"', pe.Executable, codefullpath);
    disp(cmdlinestr)
    [status, cmdout] = system(cmdlinestr, '-echo');
    if status ~= 0
        if pkg.i_isvalid(fw)
            gui.gui_waitbar(fw, true);
        end
        error('%s', cmdout);
    end
end

tmpfilelist = {'X.mat', 'g.csv', 'c.csv'};
pkg.i_deletefiles(tmpfilelist);   % always clear stale files, so a failed
% run cannot leave a previous run's output to be picked up as this one's

if issparse(sce.X)
    X = single(full(sce.X));
else
    X = single(sce.X);
end
save('X.mat','-v7.3',"X");
g = sce.g;
writetable(table(g),'g.csv','WriteVariableNames',false);
% The h5ad obs index has to be unique, and i_makeattributestable reads
% sce.c_cell_id off the object, so the ids are made unique here. SCE is a
% handle object, though, so this used to rename the caller's cell barcodes
% permanently: exporting a file silently rewrote the barcodes in the live
% dataset it was exporting. Restore them on the way out, whichever way
% this function leaves.
originalCellId = sce.c_cell_id;
restoreCellId = onCleanup(@() i_restorecellid(sce, originalCellId));
sce.c_cell_id = matlab.lang.makeUniqueStrings(sce.c_cell_id);
T = pkg.i_makeattributestable(sce);
writetable(T,'c.csv');
% disp('Files written.');

if verbose && pkg.i_isvalid(fw)
    gui.gui_waitbar(fw, [], [], 'Checking Python environment is complete');
    pause(0.5);
    gui.gui_waitbar(fw, [], [], 'Running py\_writeh5ad...');
end

codefullpath = fullfile(codepth,'script.py');
pkg.i_addwd2script(codefullpath, wkdir, 'python');
cmdlinestr = sprintf('"%s" "%s"', pe.Executable, codefullpath);
disp(cmdlinestr)
[status1] = system(cmdlinestr, '-echo');
[status2] = movefile('output.h5ad',fname);

if status1 == 0 && status2 == 1
    succeeded = true;
    if verbose && pkg.i_isvalid(fw)
        gui.gui_waitbar(fw, false, 'File is written.');
    end
else
    if verbose && pkg.i_isvalid(fw)
        gui.gui_waitbar(fw, true, 'File is failed to save.');
    end
end

if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
end

function i_restorecellid(sce, cellid)
% Put back the barcodes made unique for the h5ad obs index.
sce.c_cell_id = cellid;
end
