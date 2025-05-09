function [pred] = r_copykat(sce, wkdir)

pred = [];
if nargin < 2 
   wkdir = tempdir;
end
% PMID: 33462507


isdebug = true;
oldpth = pwd();
[isok, msg, codepath] = commoncheck_R('R_copykat');
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
end
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
cd(oldpth);
end



