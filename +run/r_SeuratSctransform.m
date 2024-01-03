function [X] = r_SeuratSctransform(X, genelist, wkdir)

if nargin < 3, wkdir = tempdir; end

isdebug = false;
oldpth = pwd();
[isok, msg, codepath] = commoncheck_R('R_SeuratSctransform');
if ~isok, error(msg);
    return;
end
if ~isempty(wkdir) && isfolder(wkdir), cd(wkdir); end

tmpfilelist = {'input.mat', 'output.h5', 'input.txt', 'output.txt', 'g.txt'};
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
lastwarn('')
save('input.mat', 'X', '-v7.3');
writematrix(genelist, 'g.txt');
[warnMsg, warnId] = lastwarn;
if ~isempty(warnMsg)
    disp(warnId)
    if exist('./input.mat', 'file'), delete('./input.mat'); end
    disp('Writing data into input.txt...')
    sc_writefile('input.txt', X, genelist);
end


Rpath = getpref('scgeatoolbox', 'rexecutablepath');
codefullpath = fullfile(codepath,'script.R');
pkg.RunRcode(codefullpath, Rpath);

if exist('output.h5', 'file')
    X = h5read('output.h5', '/X');
elseif exist('output.txt', 'file')
    X = readmatrix('output.txt');
end
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
cd(oldpth);
end
