function [y] = py_ikarus(X, g, wkdir)
if nargin < 3, wkdir = 'C:/Users/jcai/Downloads'; end

isdebug = true;
y = [];
prgfoldername = 'py_ikarus';

oldpth = pwd();
[pyok, ~, x] = run.pycommon(prgfoldername);
if ~pyok, return; end
pw1 = fileparts(mfilename('fullpath'));
codepth = fullfile(pw1, 'external', prgfoldername);
if isempty(wkdir) || ~isfolder(wkdir)
    wkdir = codepth;
    cd(codepth);
else
    disp('Using working directory provided.');
    cd(wkdir);
    copyfile(fullfile(codepth,'signatures.gmt'),'.\');
    copyfile(fullfile(codepth,'core_model.joblib'),'.\');
end

tmpfilelist = {'X.mat', 'g.csv'};
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
if issparse(X), X = full(X); end
save('X.mat', '-v7.3', 'X');
g = upper(g);
writetable(table(g),'g.csv','WriteVariableNames',false);
disp('Input files written.');

codefullpath = fullfile(codepth,'script.py');
pkg.i_addwd2script(codefullpath, wkdir, 'python');

cmdlinestr = sprintf('"%s" "%s"', x.Executable, codefullpath);
disp(cmdlinestr)
[status] = system(cmdlinestr, '-echo');

% [status] = run.pycommon2(x, codepth, prgfoldername);

if status == 0 && exist('out/test/scores.csv', 'file')
    y = readtable("out/test/scores.csv");
end

if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
cd(oldpth);
end
