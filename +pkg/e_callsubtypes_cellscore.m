function [sce]=e_callsubtypes_cellscore(sce,targettype,k)

if nargin<2, targettype="T cells"; end
if nargin<3, k=5; end
% ref: https://www.biocompare.com/Editorial-Articles/569888-A-Guide-to-T-Cell-Markers/
% ref: https://www.cellsignal.com/pathways/immune-cell-markers-human

c_cell_type_tx=erase(sce.c_cell_type_tx,"_{"+digitsPattern+"}");

if ~ismember(targettype,c_cell_type_tx)
    warning('Target cell type not found.');
    return;
end

sce2=removecells(sce,~ismember(c_cell_type_tx,targettype));

if length(unique(sce2.c_cluster_id))==1
    id=sc_cluster_s(sce2.s,k);
else
    id=sce2.c_cluster_id;
end
[id]=grp2idx(id);

% [T]=run.mt_alona(sce2.X,sce2.g,id);
    
posg1=["Cd3d","Cd3e","Cd3g","Cd8a","Cd8a1"]; 
negg1=["Gzma","Gzmb","Pdcd1","Ctla4"];
posg2=["Cd3d","Cd3e","Cd3g","Cd4","Il2ra","Foxp3"];
negg2=[];
posg3=["Cd3d","Cd3e","Cd3g","Cd8a","Cd8a1","Gzma","Gzmb","Pdcd1","Ctla4"];
negg3=[];

posg4=["Cd3g"];
negg4=["Cd4","Cd8a"];


cs1=sc_cellscore_admdl(sce2.X,sce2.g,posg1,negg1);
cs2=sc_cellscore_admdl(sce2.X,sce2.g,posg2,negg2);
cs3=sc_cellscore_admdl(sce2.X,sce2.g,posg3,negg3);
cs4=sc_cellscore_admdl(sce2.X,sce2.g,posg4,negg4);
cstype=["T\_CD8","Treg","Tex\_CD8","Double\_negative_Treg"].';
[~,idx] = max([cs1,cs2,cs3,cs4],[],2);
a=cstype(idx);
% CD8+ T_CD8 Cd3d, cd3e, cd3g, Cd8a Cd8a1
% CD4+ regulatory Treg Cd3d, cd3e, cd3g, cd4 Il2ra, Foxp3
% CD4+ T_CD4 Cd3d, cd3e, cd3g, cd4
% CD8+ exhausted Tex_CD8   Cd3d, cd3e, cd3g, Cd8a cd8b1 Gzma Gzmb Gzmk Pdcd1 Ctla4
% Double negative T_DN Cd3d cd3e Cd3g Il17a Pdcd1
[c,cL]=grp2idx(a);
cL=cL(grpstats(c,id,@mode));
annolabels=cL(id);
%annolabels=strrep(annolabels,"_","\_");
sce.c_cell_type_tx(ismember(c_cell_type_tx,targettype))=annolabels;



