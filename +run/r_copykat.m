function [pred] = r_copykat(sce, wkdir, speciesid, isdebug)

pred = [];
if nargin < 4, isdebug = true; end
if nargin < 3
    speciesid = 'human';   % mouse
end
if nargin < 2
   wkdir = pkg.i_tempdirfile();
end
% PMID: 33462507


oldpth = pwd();
cleanupCwd = onCleanup(@() cd(oldpth));
[isok, msg, codepath] = commoncheck_R('R_copykat');
if ~isok, error('%s', msg); end
if ~isempty(wkdir) && isfolder(wkdir), cd(wkdir); end

tmpfilelist = {'input.h5'};
pkg.i_deletefiles(tmpfilelist);   % always clear stale files, so a failed
% run cannot leave a previous run's output to be picked up as this one's
pkg.e_writeh5(full(sce.X), sce.g, 'input.h5');
Rpath = getpref('scgeatoolbox', 'rexecutablepath',[]);
if isempty(Rpath)
    error('R environment has not been set up.');
end

if string(lower(speciesid))=="human"
    codefullpath = fullfile(codepath,'script_hg20.R');
elseif string(lower(speciesid))=="mouse"
    codefullpath = fullfile(codepath,'script_mm10.R');
else
    codefullpath = fullfile(codepath,'script.R');
end
pkg.i_runrcode(codefullpath, Rpath);
pause(3);
outfile = "test_copykat_prediction.txt";
if exist(outfile,'file')
    t = readtable(outfile, ...
        "ReadVariableNames", true, "Delimiter",'\t',...
        "VariableNamingRule", "modify");
    assert(height(t)==sce.NumCells);
    idx = str2double(extractAfter(string(t.cell_names),1));
    [~, sortid] = sort(idx);
    pred = string(t.copykat_pred);
    pred = pred(sortid);
    % y = zeros(sce.NumCells, 1);
    % y(pred == "aneuploid") = 1;
    % y(pred == "diploid") = 0;
    % y = y(sortid);
else
    error('run.r_copykat:noOutput', ...
        ['R finished but did not write %s to %s. The R console ', ...
         'output above should say why.'], outfile, pwd);
end
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
end
