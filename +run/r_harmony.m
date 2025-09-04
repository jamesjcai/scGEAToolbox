function [sout] = r_harmony(s, batchid, wkdir, isdebug)
arguments
    s(:, :) {mustBeNumeric}
    batchid
    wkdir = []
    isdebug = false
end

if isempty(wkdir), wkdir = tempdir; end

oldpth = pwd();
[isok, msg, codepth] = commoncheck_R('R_harmony');
if ~isok, error(msg);
    return;
end
if ~isempty(wkdir) && isfolder(wkdir), cd(wkdir); end


sout = [];
tmpfilelist = {'input.h5', 'output.h5'};
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end

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
pkg.RunRcode(codefullpath, Rpath);

if exist('output.h5', 'file')
    % load("output.mat", "sout")
    sout = h5read("output.h5", "/harmony_embeddings");
end


if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
cd(oldpth);

end
