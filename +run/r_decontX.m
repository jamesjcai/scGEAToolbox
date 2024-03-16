function [X, contamination] = r_decontX(X, wkdir, isdebug)
%Run decontX decontamination
%
% see also: run.r_SoupX
% https://cran.r-project.org/web/packages/SoupX/vignettes/pbmcTutorial.html

if nargin < 2, wkdir = tempdir; end
if nargin < 3, isdebug = false; end

oldpth = pwd();
[isok, msg, codepath] = commoncheck_R('R_decontX');

if ~isok, error(msg); end
if ~isempty(wkdir) && isfolder(wkdir), cd(wkdir); end

if isa(X, 'SingleCellExperiment')
    X = X.X;
end

tmpfilelist = {'input.mat', 'output.h5'};
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end

if issparse(X), X = full(X); end
save('input.mat', 'X', '-v7.3');
Rpath = getpref('scgeatoolbox', 'rexecutablepath',[]);
if isempty(Rpath)
    error('R environment has not been set up.');
end

codefullpath = fullfile(codepath,'script.R');
pkg.RunRcode(codefullpath, Rpath);
if exist('output.h5', 'file')
    X = h5read('output.h5', '/X');
    contamination = h5read('output.h5', '/contamination');
    % load('output.mat','X','contamination')
end

if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
cd(oldpth);
end
