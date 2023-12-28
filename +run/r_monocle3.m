function [t] = r_monocle3(X, idx, wkdir, isdebug)
% Run Monocle pseudotime analysis
% [t]=run.r_monocle3(X,idx);
if nargin < 2, idx = [1, 2]; end
if nargin < 3, wkdir = ''; end
if nargin < 4, isdebug = false; end

t = [];
oldpth = pwd();
[isok, msg, codepath] = commoncheck_R('R_monocle3');
if ~isok, error(msg); end
if ~isempty(wkdir), cd(wkdir); end

tmpfilelist = {'input.h5', 'output.h5', 'input.mat'};
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end

if issparse(X), X = full(X); end
save('input.mat', 'X', 'idx', '-v7.3');
%pkg.e_writeh5(X,[],'input.h5');
Rpath = getpref('scgeatoolbox', 'rexecutablepath');

codefullpath = fullfile(codepath,'script.R');
pkg.RunRcode(codefullpath, Rpath);
if exist('output.h5', 'file')
    t = h5read('output.h5', '/t');
end
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
cd(oldpth);
end
