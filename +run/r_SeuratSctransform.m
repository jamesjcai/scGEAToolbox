function [X,scale_X] = r_SeuratSctransform(X, genelist, wkdir)

if nargin < 3, wkdir = tempdir; end
if nargin < 2, genelist = string(1:size(X, 1)); end
isdebug = false;
oldpth = pwd();
[isok, msg, codepath] = commoncheck_R('R_SeuratSctransform');
if ~isok, error(msg);
    return;
end
if ~isempty(wkdir) && isfolder(wkdir), cd(wkdir); end

tmpfilelist = {'input.mat', 'output.h5', 'input.txt', 'output.txt', ...
    'g.txt', 'output_data.txt', 'output_scale_data.txt'};
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
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
pkg.RunRcode(codefullpath, Rpath);

if exist('output.h5', 'file')
    X = h5read('output.h5', '/data');
elseif exist('output_data.txt', 'file')
    X = readmatrix('output_data.txt');
end
if nargout>1
    if exist('output.h5', 'file')
        scale_X = h5read('output.h5', '/scale_data');
    elseif exist('output_scale_data.txt', 'file')
        scale_X = readmatrix('output_scale_data.txt');
    end
end
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
cd(oldpth);
end
