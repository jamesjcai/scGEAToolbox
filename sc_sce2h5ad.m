function [succeeded] = sc_sce2h5ad(sce, filename, verbose, pe)
% Write SCE to H5AD file
%
% see also: sc_readh5adfile

if nargin < 4, pe = []; end
if nargin < 3, verbose = true; end

succeeded = false;
if nargin < 2
    [filename, pathname] = uiputfile({'*.h5ad'; '*.*'}, 'Save as');
    if ~(filename), return; end
    filename = fullfile(pathname, filename);
end

wkdir = pkg.i_tempdirfile();
try
    [succeeded] = run.py_writeh5ad(sce, filename, wkdir, false, verbose, pe);
catch ME
    errordlg(ME.message);
end
end
