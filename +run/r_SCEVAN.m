function [s, t] = r_SCEVAN(sce, wkdir, ~, species)

if nargin<3, SUBCLONES = false; end
if nargin<4, species = 'human'; end
s = []; t = [];
if nargin < 2 
   wkdir = tempdir;
end
% PMID: 36841879


isdebug = true;
oldpth = pwd();
[isok, msg, codepath] = commoncheck_R('R_SCEVAN');
if ~isok, error(msg);
    return;
end
if ~isempty(wkdir) && isfolder(wkdir), cd(wkdir); end

tmpfilelist = {'input.h5','output.csv'};
pkg.i_deletefiles(tmpfilelist);
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
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
pkg.RunRcode(codefullpath, Rpath);
pause(3);
outfile = "output.csv";
if exist(outfile,'file')
    t = readtable(outfile);
    assert(height(t)==sce.NumCells);
    s = string(t.class);
end
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
cd(oldpth);
end

