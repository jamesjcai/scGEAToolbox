function py_writeh5ad(sce, fname)
isdebug=fales;

prgfoldername = 'py_writeh5ad';
[pyok, wrkpth, x] = run.pycommon(prgfoldername);
if ~pyok, return; end
tmpfilelist = {'X.h5', 'obs.csv', 'var.csv'};
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end


[status] = run.pycommon2(x, wrkpth, prgfoldername);

oldpth = pwd();
pw1 = fileparts(mfilename('fullpath'));
wrkpth = fullfile(pw1, 'external', 'py_writeh5ad');
cd(wrkpth);
tmpfilelist = {'X.h5', 'obs.csv', 'var.csv'};
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end

[status] = run.pycommon2(x, wrkpth, prgfoldername);

if status == 0 && exist('output.txt', 'file')

end

if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
cd(oldpth);

end
