function py_datamapplot(sce, thisc, wkdir, isdebug)

if nargin < 4, isdebug = true; end
if nargin < 3, wkdir = []; end
if nargin < 2, thisc = sce.c_cell_type_tx; end

oldpth = pwd();
pw1 = fileparts(mfilename('fullpath'));
codepth = fullfile(pw1, 'external', 'py_datamapplot');

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
%cmdlinestr = sprintf('"%s" "%s%srequire.py"', ...
%    x.Executable, codepth, filesep);
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

tmpfilelist = {'input.mat', 'c.txt'};
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end

writematrix(thisc, 'c.txt', 'Delimiter', ',');
s = sce.s;
if issparse(s), s = full(s); end
save('input.mat', '-v7.3', 's');
disp('Input files written.');

if isvalid(fw)
    gui.gui_waitbar(fw, [], [], 'Checking Python environment is complete');
    pause(0.5);
    gui.gui_waitbar(fw, [], [], 'Running DataMapPlot...');
end
% fw = gui.gui_waitbar([],[],'Running DataMapPlot...');
codefullpath = fullfile(codepth,'script.py');

pkg.i_addwd2script(codefullpath, wkdir, 'python');

cmdlinestr = sprintf('"%s" "%s"', x.Executable, codefullpath);
disp(cmdlinestr)
[status] = system(cmdlinestr, '-echo');

if status == 0 && isvalid(fw)
    gui.gui_waitbar(fw, [], 'DataMapPlot is complete');
end

if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
cd(oldpth);
end
