function [sce]=e_callsubtypes_alona(sce,targettype,k)

if nargin<2, targettype="T cells"; end
if nargin<3, k=[]; end
% ref: https://www.biocompare.com/Editorial-Articles/569888-A-Guide-to-T-Cell-Markers/
% ref: https://www.cellsignal.com/pathways/immune-cell-markers-human

c_cell_type_tx=erase(sce.c_cell_type_tx,"_{"+digitsPattern+"}");

if ~ismember(targettype,c_cell_type_tx)
    warning('Target cell type not found.');
    return;
end

annolabels=sce.c_cell_type_tx(ismember(c_cell_type_tx,targettype));

sce2=removecells(sce,~ismember(c_cell_type_tx,targettype));

if length(unique(sce2.c_cluster_id))==1 || ~isempty(k)
    id=sc_cluster_s(sce2.s,k);
else
    id=sce2.c_cluster_id;
end
[clusterid]=grp2idx(id);

switch targettype
    case 'T cells'
        targettag='tcells';
    case 'Neurons'
        targettag='neurons';
    case 'Fibroblasts'
        targettag='fibroblasts';
end

for ix=1:max(clusterid)
    [T]=run.alona(sce2.X(:,clusterid==ix),sce2.g,[],...
        'subtype',targettag);
    ctxt = sprintf('%s_{%d}', T.C1_Cell_Type{1}, ix);
    annolabels(clusterid==ix)=ctxt;
end
% [c,cL]=grp2idx(a);
% cL=cL(grpstats(c,id,@mode));
% annolabels=cL(id);
% %annolabels=strrep(annolabels,"_","\_");
sce.c_cell_type_tx(ismember(c_cell_type_tx,targettype))=annolabels;



