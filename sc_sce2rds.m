function [status]=sc_sce2rds(sce,filename)
%Write SCE to Seurat/RDS file
%
%see also: sc_readrdsfile

status=0;
if nargin < 2
    [filename, pathname] = uiputfile( ...
       {'*.*'},'Save as');
	if ~(filename), return; end
	filename=[pathname,filename];
end
status=run.r_saveSeuratRds(sce,filename);
