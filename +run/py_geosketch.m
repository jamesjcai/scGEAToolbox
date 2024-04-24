function [idx] = py_geosketch(X, n)

isdebug = false;
if nargin < 2
    n = min([2000, round(0.75*size(X, 2))]);
end
idx = [];

prgfoldername = 'py_geosketch';

oldpth = pwd();
[pyok, wrkpth, x] = run.pycommon(prgfoldername);
if ~pyok, return; end

tmpfilelist = {'input.mat', 'output.mat'};
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
if issparse(X), X = full(X); end
save('input.mat', '-v7.3', 'X', 'n');
disp('Input file written.');


%     fw=gui.gui_waitbar([],[],'Running geosketch...');
%     cmdlinestr=sprintf('"%s" "%s%sscript.py"', ...
%         x.Executable,wrkpth,filesep);
%     disp(cmdlinestr)
%     [status]=system(cmdlinestr,'-echo');
%     if isvalid(fw)
%         gui.gui_waitbar(fw,[],'Running geosketch is complete');
%     end

[status] = run.pycommon2(x, wrkpth, prgfoldername);

if status == 0 && exist('output.mat', 'file')
    load("output.mat", "idx")
    idx = idx + 1;
end

if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
cd(oldpth);
end
