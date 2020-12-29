function [sce]=sc_mergensces(sces,method)

if nargin<2, method='intersect'; end
if ~iscell(sces), error('SCES is not a cell array'); end
sce=sces{1};
c=ones(sce.NumCells,1);
for k=2:length(sces)
    c=[c; k*ones(sces{k}.NumCells,1)];
    [sce]=sc_merge2sces(sce,sces{k},method);
end
sce.c_batch_id=c;
end
