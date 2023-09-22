function [T] = py_GenKI(X, g, idx)

prgfoldername = 'py_GenKI';
narginchk(3, 3);
assert(size(X, 1) == length(g));

T = [];
isdebug = false;
oldpth = pwd();
[pyok, wrkpth, x] = run.pycommon(prgfoldername);
if ~pyok
    error('GenKI (requires Python) has not been installed or set up correctly.')
end

%     pw1=fileparts(mfilename('fullpath'));
%     wrkpth=fullfile(pw1,'external','py_GenKI');
%     cd(wrkpth);
%
%     fw = gui.gui_waitbar([],[],'Checking Python environment...');
%     x=pyenv;
%     try
%         pkg.i_add_conda_python_path;
%     catch
%
%     end
%     cmdlinestr=sprintf('"%s" "%s%srequire.py"', ...
%             x.Executable,wrkpth,filesep);
%     disp(cmdlinestr)
%     [status,cmdout]=system(cmdlinestr,'-echo');
%     if status~=0
%         gui.gui_waitbar(fw,true);
%         cd(oldpth);
%         waitfor(errordlg(sprintf('%s',cmdout)));
%         error('Python GenKI has not been installed properly.');
%     end
%
%     if isvalid(fw)
%         gui.gui_waitbar(fw,[],'Checking Python environment is complete');
%     end

tmpfilelist = {'X.mat', 'g.txt', 'c.txt', ...
    'idx.mat', 'output.csv', fullfile('GRNs', 'pcNet_example.npz')};
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end

save('X.mat', '-v7.3', 'X');
save('idx.mat', '-v7.3', 'idx');
writematrix(g, 'g.txt');
writematrix(ones(size(X, 2), 1), 'c.txt');
disp('Input X g c written.');


% fw=gui.gui_waitbar([],[],'Running GenKI...');

%    cmdlinestr=sprintf('"%s" "%s%sscript.py"', ...
%        x.Executable,wrkpth,filesep);
%    disp(cmdlinestr)
%    [status]=system(cmdlinestr,'-echo');

[status] = run.pycommon2(x, wrkpth, prgfoldername);

% if isvalid(fw)
%     gui.gui_waitbar(fw,[],'Running GenKI is complete');
% end

if status == 0 && exist('output.csv', 'file')
    T = readtable('output.csv');
    T.Properties.VariableNames{1} = 'gene';
end

if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
cd(oldpth);
end
