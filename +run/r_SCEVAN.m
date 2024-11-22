function r_SCEVAN(sce, wkdir)

% PMID: 36841879

if nargin < 2, wkdir = tempdir; end

isdebug = true;
oldpth = pwd();
[isok, msg, codepath] = commoncheck_R('R_SCEVAN');
if ~isok, error(msg);
    return;
end
if ~isempty(wkdir) && isfolder(wkdir), cd(wkdir); end

tmpfilelist = {'input.h5'};
pkg.i_deletefiles(tmpfilelist);
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
pkg.e_writeh5(full(sce.X), sce.g, 'input.h5');
Rpath = getpref('scgeatoolbox', 'rexecutablepath',[]);
if isempty(Rpath)
    error('R environment has not been set up.');
end
codefullpath = fullfile(codepath,'script.R');
pkg.RunRcode(codefullpath, Rpath);
pause(3);
outfile = "output.txt";
if exist(outfile,'file')
    warning off
    t = readtable(outfile, ...
        "ReadVariableNames", true, "Delimiter",'\t',...
        "VariableNamingRule", "modify");
    warning on
    assert(height(t)==sce.NumCells);
end
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
cd(oldpth);
end



