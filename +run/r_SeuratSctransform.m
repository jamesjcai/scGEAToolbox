function [X,scale_X] = r_SeuratSctransform(X, genelist, wkdir)

if nargin < 3, wkdir = pkg.i_tempdirfile(); end
if nargin < 2, genelist = string(1:size(X, 1)); end
isdebug = false;
oldpth = pwd();
cleanupCwd = onCleanup(@() cd(oldpth));
[isok, msg, codepath] = commoncheck_R('R_SeuratSctransform');
if ~isok, error('%s', msg); end
if ~isempty(wkdir) && isfolder(wkdir), cd(wkdir); end

tmpfilelist = {'input.mat', 'output.h5', 'input.txt', 'output.txt', ...
'g.txt', 'output_data.txt', 'output_scale_data.txt'};
pkg.i_deletefiles(tmpfilelist);   % always clear stale files, so a failed
% run cannot leave a previous run's output to be picked up as this one's
lastwarn('')
if issparse(X), X = full(X); end
save('input.mat', 'X', '-v7.3');
writematrix(genelist, 'g.txt');
[warnMsg, warnId] = lastwarn;
if ~isempty(warnMsg)
    disp(warnId)
    if exist('./input.mat', 'file'), delete('./input.mat'); end
    disp('Writing data into input.txt...')
    sc_writefile('input.txt', X, genelist);
end


Rpath = getpref('scgeatoolbox', 'rexecutablepath',[]);
if isempty(Rpath)
    error('R environment has not been set up.');
end
codefullpath = fullfile(codepath,'script.R');
pkg.i_addwd2script(codefullpath, wkdir, 'R');
pkg.i_runrcode(codefullpath, Rpath);

if exist('output.h5', 'file')
    X = h5read('output.h5', '/data');
elseif exist('output_data.txt', 'file')
    X = readmatrix('output_data.txt');
else
    error('run.r_SeuratSctransform:noOutput', ...
        ['R finished but wrote neither output.h5 nor output_data.txt to ', ...
         '%s. The R console output above should say why.'], pwd);
end
if nargout>1
    if exist('output.h5', 'file')
        scale_X = h5read('output.h5', '/scale_data');
    elseif exist('output_scale_data.txt', 'file')
        scale_X = readmatrix('output_scale_data.txt');
    else
        error('run.r_SeuratSctransform:noScaleOutput', ...
            ['R produced the normalised matrix but no scaled matrix in ', ...
             '%s. Call with one output if the scaled matrix is not needed.'], pwd);
    end
end
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
end
