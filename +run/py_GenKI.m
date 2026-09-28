function [T] = py_GenKI(X, g, idx, wkdir, isdebug)

if nargin < 5, isdebug = true; end
if nargin < 4, wkdir = pkg.i_tempdirfile(); end

oldpth = pwd();
cleanupCwd = onCleanup(@() cd(oldpth));
pw1 = fileparts(mfilename('fullpath'));
codepth = fullfile(pw1, '..', 'external', 'py_GenKI');

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

codepth = pkg.i_normalizepath(codepth);

codefullpath = fullfile(codepth,'require.py');
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


try
    tmpfilelist = {'X.mat', 'g.txt', 'c.txt', 'pcnet_Source.mat', ...
        'idx.mat', 'output.csv', fullfile('GRNs', 'pcNet_example.npz')};
    pkg.i_deletefiles(tmpfilelist);   % always clear stale files, so a failed
    % run cannot leave a previous run's output to be picked up as this one's
    if issparse(X)
        X = single(full(X));
    else
        X = single(X);
    end
    save('X.mat', '-v7.3', 'X');
    save('idx.mat', '-v7.3', 'idx');
    writematrix(g, 'g.txt');
    writematrix(ones(size(X, 2), 1), 'c.txt');
catch ME
    if pkg.i_isvalid(fw)
         gui.myWaitbar([], fw, true);
    end
    errordlg(ME.message,'');
    return;
end
if pkg.i_isvalid(fw)
    gui.myWaitbar([], fw, [], [], 'Checking Python environment is complete');
    pause(0.5);
    gui.myWaitbar([], fw, [], [], 'Running GenKI...');
end


if pkg.i_isvalid(fw)
    gui.myWaitbar([], fw, [], [], 'Building pcnet_Source network...');
end
A1 = net.pcrnet(X, 3, false, true, false, false, pkg.i_usegpu(X));
A1 = A1 ./ max(abs(A1(:)));
A = ten.e_filtadjc(A1, 0.75, false);
save('pcnet_Source.mat', 'A', '-v7.3');
if pkg.i_isvalid(fw)
    gui.myWaitbar([], fw, [], [], 'pcnet_Source.mat saved.');
end

codefullpath = fullfile(codepth,'script.py');
pkg.i_addwd2script(codefullpath, wkdir, 'python');
cmdlinestr = sprintf('"%s" "%s"', x.Executable, codefullpath);
disp(cmdlinestr)
[status] = system(cmdlinestr, '-echo');

if status == 0 && pkg.i_isvalid(fw)
    gui.myWaitbar([], fw, [], 'py_GenKI is complete');
end


%    cmdlinestr=sprintf('"%s" "%s%sscript.py"', ...
%        x.Executable,wrkpth,filesep);
%    disp(cmdlinestr)
%    [status]=system(cmdlinestr,'-echo');

if status == 0 && exist('output.csv', 'file')
    T = readtable('output.csv');
    T.Properties.VariableNames{1} = 'gene';
end

if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
end
