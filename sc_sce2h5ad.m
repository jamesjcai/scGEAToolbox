function [status] = sc_sce2h5ad(sce, filename)
%Write SCE to H5AD file
%
%see also: sc_readh5adfile

status = 1;
if nargin < 2
     [filename, pathname] = uiputfile({'*.h5ad'; '*.*'}, 'Save as');
    if ~(filename), return; end
    filename = [pathname, filename];
end

%extprogname = 'py_writeh5ad';
%preftagname = 'externalwrkpath';
%[wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname);
wkdir = tempdir;
try
    run.py_writeh5ad(sce, filename, wkdir, false);
catch ME
    errordlg(ME.message);
end
status = 0;
