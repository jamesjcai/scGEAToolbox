function [isDoublet, doubletscore] = py_scrublet_new(X, wkdir, isdebug)

if nargin < 3, isdebug = true; end
if nargin < 2, wkdir = []; end

isDoublet = [];
doubletscore = [];
prgfoldername = 'py_scrublet';
% [pyok, wrkpth, x] = run.pycommon(prgfoldername);
% if ~pyok, return; end
% isdebug = false;

oldpth = pwd();
pw1 = fileparts(mfilename('fullpath'));
codepth = fullfile(pw1, 'external', prgfoldername);

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





tmpfilelist = {'input.mat', 'output.mat'};
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end

if issparse(X), X = full(X); end

save('input.mat', '-v7.3', 'X');
disp('Input file written.');

if isvalid(fw)
    gui.gui_waitbar(fw, [], [], 'Checking Python environment is complete');
    pause(0.5);
    gui.gui_waitbar(fw, [], [], 'Running Scrublet...');
end
% fw = gui.gui_waitbar([],[],'Running Scrublet...');
codefullpath = fullfile(codepth,'script.py');

pkg.i_addwd2script(codefullpath, wkdir, 'python');

cmdlinestr = sprintf('"%s" "%s"', x.Executable, codefullpath);
disp(cmdlinestr)
[status] = system(cmdlinestr, '-echo');



% cmdlinestr=sprintf('"%s" "%s%sscript.py"', ...
%     x.Executable,wrkpth,filesep);
% disp(cmdlinestr)
% [status]=system(cmdlinestr,'-echo');
% [status] = run.pycommon2(x, wrkpth, prgfoldername);

if status == 0 && exist('output.mat', 'file')
    load("output.mat", 'isDoublet', 'doubletscore')
end
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
cd(oldpth);

gui.gui_waitbar(fw);
end
