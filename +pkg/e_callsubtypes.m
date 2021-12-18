function [sce2]=e_callsubtypes(sce,targettype)
if nargin<2, targettype="T cells"; end
c_cell_type_tx=erase(sce.c_cell_type_tx,"_{"+digitsPattern+"}");
sce2=[];
if ~ismember(targettype,c_cell_type_tx)
    warning('Target cell type not found.');
    return;
end

sce2=removecells(sce,~ismember(c_cell_type_tx,targettype));

if length(unique(sce2.c_cluster_id))==1
    id=sc_cluster_s(sce2.s,5);
else
    id=sce2.c_cluster_id;
end
    [T]=run.alona(sce2.X,sce2.g,id);



% CD8+ T_CD8 Cd3d, cd3e, cd3g, Cd8a Cd8a1
% CD4+ regulatory Treg Cd3d, cd3e, cd3g, cd4 Il2ra, Foxp3
% CD4+ T_CD4 Cd3d, cd3e, cd3g, cd4
% CD8+ exhausted Tex_CD8   Cd3d, cd3e, cd3g, Cd8a cd8b1 Gzma Gzmb Gzmk Pdcd1 Ctla4
% Double negative T_DN Cd3d cd3e Cd3g Il17a Pdcd1

