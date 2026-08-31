function [s, t] = r_SCEVAN(sce, wkdir, ~, species, isdebug)

if nargin<3, SUBCLONES = false; end
if nargin<4, species = 'human'; end
if nargin<5, isdebug = true; end
s = []; t = [];
if nargin < 2
   wkdir = pkg.i_tempdirfile();
end
% PMID: 36841879


oldpth = pwd();
cleanupCwd = onCleanup(@() cd(oldpth));
[isok, msg, codepath] = commoncheck_R('R_SCEVAN');
if ~isok, error('%s', msg); end
if ~isempty(wkdir) && isfolder(wkdir), cd(wkdir); end

tmpfilelist = {'input.h5','output.csv'};
pkg.i_deletefiles(tmpfilelist);
pkg.i_deletefiles(tmpfilelist);   % always clear stale files, so a failed
% run cannot leave a previous run's output to be picked up as this one's
pkg.e_writeh5(full(sce.X), sce.g, 'input.h5');
Rpath = getpref('scgeatoolbox', 'rexecutablepath',[]);
if isempty(Rpath)
    error('R environment has not been set up.');
end
if strcmpi(species, 'human')
    codefullpath = fullfile(codepath,'script_human.R');
else
    codefullpath = fullfile(codepath,'script_mouse.R');
end
pkg.i_runrcode(codefullpath, Rpath);
pause(3);
outfile = "output.csv";
if exist(outfile,'file')
    t = readtable(outfile);
    assert(height(t)==sce.NumCells);
    s = string(t.class);
else
    error('run.r_SCEVAN:noOutput', ...
        ['R finished but did not write %s to %s. The R console ', ...
         'output above should say why.'], outfile, pwd);
end
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
end
