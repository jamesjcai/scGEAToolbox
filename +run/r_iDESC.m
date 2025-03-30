function [T] = r_iDESC(sce, group_var, subject_var, wkdir, isdebug)

% https://github.com/yl883/iDESC/

% \item{subject_var}{The name of subject information in meta}
% \item{group_var}{The name of group/disease information for DE analysis in meta}

if nargin < 4, wkdir = tempdir; end
if nargin < 5, isdebug = false; end

T = [];
oldpth = pwd();
[isok, msg, codepth] = commoncheck_R('R_iDESC');
if ~isok, error(msg); end
if ~isempty(wkdir) && isfolder(wkdir), cd(wkdir); end

tmpfilelist = {'g.txt', 'a.txt', 'input.mat'};
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end

if issparse(X), X = full(X); end
save('input.mat', 'X', '-v7.3');

writematrix(sce.g, 'g.txt');
t = pkg.makeattributestable(sce);
sequencing_depth = sum(X)';
t = [t table(sequencing_depth, group_var, subject_var)];
writetable(t,'a.txt');

Rpath = getpref('scgeatoolbox', 'rexecutablepath',[]);
if isempty(Rpath)
    error('R environment has not been set up.');
end
codefullpath = fullfile(codepth,'script.R');
pkg.i_addwd2script(codefullpath, wkdir, 'R');
pkg.RunRcode(codefullpath, Rpath);
if exist('output.csv', 'file')
    T = readtable('output.csv');
end
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
cd(oldpth);
end
