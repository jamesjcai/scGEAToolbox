function [sout] = r_harmony(s, batchid, wkdir, isdebug)
arguments
s(:, :) {mustBeNumeric}
batchid
wkdir = []
isdebug = false
end

if isempty(wkdir), wkdir = pkg.i_tempdirfile(); end

oldpth = pwd();
cleanupObj = onCleanup(@() cd(oldpth));
[isok, msg, codepth] = commoncheck_R('R_harmony');
if ~isok, error('%s', msg); end
if ~isempty(wkdir) && isfolder(wkdir), cd(wkdir); end


sout = [];
tmpfilelist = {'input.h5', 'output.h5'};
pkg.i_deletefiles(tmpfilelist);   % always clear stale files, so a failed
% run cannot leave a previous run's output to be picked up as this one's

if issparse(s), s = full(s); end
filename = "input.h5";

h5create(filename, "/s", size(s));
h5write(filename, "/s", s);
batchid = string(batchid);
h5create(filename, '/batchid', size(batchid), 'Datatype', 'string');
h5write(filename, '/batchid', batchid);


Rpath = getpref('scgeatoolbox', 'rexecutablepath',[]);
if isempty(Rpath)
    error('R environment has not been set up.');
end

codefullpath = fullfile(codepth,'script.R');
pkg.i_addwd2script(codefullpath, wkdir, 'R');
pkg.i_runrcode(codefullpath, Rpath);

if exist('output.h5', 'file')
    % load("output.mat", "sout")
    sout = h5read("output.h5", "/harmony_embeddings");
else
    error('run.r_harmony:noOutput', ...
        ['R finished but did not write %s to %s. The R console ', ...
         'output above should say why.'], 'output.h5', pwd);
end


if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
% cd(oldpth);

end
