function [OUT1,OUT2,Tok]=talklr2(sce)

% https://doi.org/10.1101/2020.02.01.930602

if ~isa(sce, 'SingleCellExperiment'), error(''); end
if numel(unique(sce.c_batch_id))~=2, error(''); end
b=grp2idx(sce.c_batch_id);

[M1,M2,OUT1,OUT2]=run.i_talkr2(sce.X,sce.g,sce.c_cell_type_tx,b);
M=M2.*log2(M2./M1);
M(isnan(M))=0;
KL=real(sum(M,2));
Tok=[OUT1.Tok, table(KL)];
[Tok,idx]=sortrows(Tok,'KL','descend');

OUT1.ligand_mat=OUT1.ligand_mat(idx,:);
OUT1.receptor_mat=OUT1.receptor_mat(idx,:);
OUT1.ligandok=OUT1.ligandok(idx);
OUT1.receptorok=OUT1.receptorok(idx);
OUT1.Tok=OUT1.Tok(idx,:);

OUT2.ligand_mat=OUT2.ligand_mat(idx,:);
OUT2.receptor_mat=OUT2.receptor_mat(idx,:);
OUT2.ligandok=OUT2.ligandok(idx);
OUT2.receptorok=OUT2.receptorok(idx);
OUT2.Tok=OUT2.Tok(idx,:);

%%
for k=1:min([5 size(OUT1.ligand_mat,1)])    
    gui.i_crosstalkgraph2(OUT1,OUT2,k)
end
end