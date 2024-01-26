function [t, s, m, q] = r_monocle3(X, idx, ndim, wkdir, isdebug)
% Run Monocle pseudotime analysis
% [t]=run.r_monocle3(X,idx);
if nargin < 2, idx = [1, 2]; end
if nargin < 3 || isempty(ndim), ndim = 2; end
if nargin < 4, wkdir = tempdir; end
if nargin < 5, isdebug = false; end

t = []; s = []; m = [];
oldpth = pwd();
[isok, msg, codepth] = commoncheck_R('R_monocle3');
if ~isok, error(msg); end
if ~isempty(wkdir) && isfolder(wkdir), cd(wkdir); end

tmpfilelist = {'input.h5', 'output.h5', 'input.mat'};
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end

if issparse(X), X = full(X); end
save('input.mat', 'X', 'idx', 'ndim', '-v7.3');
%pkg.e_writeh5(X,[],'input.h5');
Rpath = getpref('scgeatoolbox', 'rexecutablepath');

codefullpath = fullfile(codepth,'script.R');
pkg.RunRcode(codefullpath, Rpath);
if exist('output.h5', 'file')
    t = h5read('output.h5', '/t');
    s = h5read('output.h5', '/s');
    m = h5read('output.h5', '/m');
    q = h5read('output.h5', '/q');
end
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
cd(oldpth);
end
