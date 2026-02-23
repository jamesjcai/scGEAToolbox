function [X, contamination] = r_decontX(X, wkdir, isdebug)
%Run decontX decontamination
%
% see also: run.r_SoupX
% https://cran.r-project.org/web/packages/SoupX/vignettes/pbmcTutorial.html

if nargin < 2, wkdir = tempdir; end
if nargin < 3, isdebug = true; end

oldpth = pwd();
cleanupObj = onCleanup(@() cd(oldpth));
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
pkg.i_addwd2script(codefullpath, wkdir, 'R');
pkg.RunRcode(codefullpath, Rpath);

outputfile1 = fullfile(wkdir,'output.h5');
outputfile2 = fullfile(codepath,'output.h5');

if exist(outputfile1, 'file')
   outputfile = outputfile1;
else
    if exist(outputfile2, 'file')
        outputfile = outputfile2;
    else
        error('output.h5 missing.');
    end
end

    X = h5read(outputfile, '/X');
    contamination = h5read(outputfile, '/contamination');

if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
end
