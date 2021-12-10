function [sce]=sc_readrdsfile(filename)
%Read Seurat/RDS file
%
%see also: sc_sce2rds
sce=[];
if nargin<1
    [fname, pathname] = uigetfile( ...
                          {'*.rds', 'RDS Format Files (*.rds)'
                           '*.*',  'All Files (*.*)'}, ...
                          'Pick a Seurat RDS file');
    if ~(fname), return; end
    filename = fullfile(pathname, fname);    
end

if ~exist(filename, 'file'), return; end
[sce]=run.readSeuratRds(filename);

end
