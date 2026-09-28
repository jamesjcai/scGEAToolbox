function [sout] = py_harmonypy(s, batchid, wkdir, isdebug)
arguments
s(:, :) {mustBeNumeric}
batchid(1, :) {mustBePositive, mustBeInteger}
wkdir = pkg.i_tempdirfile()
isdebug = true
end

oldpth = pwd();
cleanupCwd = onCleanup(@() cd(oldpth));
pw1 = fileparts(mfilename('fullpath'));
codepth = fullfile(pw1, '..', 'external', 'py_harmonypy');
codepth = pkg.i_normalizepath(codepth);

if isempty(wkdir) || ~isfolder(wkdir)
    cd(codepth);
else
    disp('Using working directory provided.');
    cd(wkdir);
end

fw = gui.myWaitbar([], [], [], 'Checking Python environment...');

x = pyenv;
try
    pkg.i_add_conda_python_path;
catch
    % best-effort: fall back to default pyenv if conda path not found
end

codefullpath = fullfile(codepth, 'require.py');
% cmdlinestr = sprintf('"%s" "%s%srequire.py"', ...
%    x.Executable, codepth, filesep);
cmdlinestr = sprintf('"%s" "%s"', x.Executable, codefullpath);

disp(cmdlinestr)
[status, cmdout] = system(cmdlinestr, '-echo');
if status ~= 0
    if pkg.i_isvalid(fw)
        gui.myWaitbar([], fw, true);
    end
    error('%s', cmdout);
end

sout = [];
% isdebug = false;
%
% %if nargin<3, usepylib=false; end
% %if nargin<2, error('[s]=run.harmonypy(s,batchid)'); end
%
% oldpth = pwd();
% [pyok, wrkpth, x] = run.pycommon(prgfoldername);
% if ~pyok, return; end

tmpfilelist = {'input.mat', 'output.mat'};
%  if ~isdebug, pkg.i_deletefiles(tmpfilelist); end

pkg.i_deletefiles(tmpfilelist);   % always clear stale files, so a failed
% run cannot leave a previous run's output to be picked up as this one's

if issparse(s), s = full(s); end
save('input.mat', '-v7.3', 's', 'batchid');
disp('Input file written.');


if pkg.i_isvalid(fw)
    gui.myWaitbar([], fw, [], [], 'Checking Python environment is complete');
    gui.myWaitbar([], fw, [], [], 'Running Harmonypy...');
end
codefullpath = fullfile(codepth,'script.py');

pkg.i_addwd2script(codefullpath, wkdir, 'python');

cmdlinestr = sprintf('"%s" "%s"', x.Executable, codefullpath);
disp(cmdlinestr)
[status] = system(cmdlinestr, '-echo');


if status == 0 && exist('output.mat', 'file')
    load("output.mat", "sout")
end

if status == 0 && pkg.i_isvalid(fw)
    gui.myWaitbar([], fw, [], 'Harmonypy is complete');
end

if ~isdebug, pkg.i_deletefiles(tmpfilelist); end

end
