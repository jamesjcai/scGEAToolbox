function [status]=sc_sce2rds(sce,filename)

status=0;
if nargin < 2
    [filename, pathname] = uiputfile( ...
       {'*.*'},'Save as');
	if ~(filename), return; end
	filename=[pathname,filename];
end
status=run.saveSeuratRds(sce,filename);
