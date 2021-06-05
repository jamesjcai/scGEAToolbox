function [OUT]=talklr(X,g,c)

% https://doi.org/10.1101/2020.02.01.930602

if nargin<2, g=[]; c=[]; end
if isa(X, 'SingleCellExperiment')    
    g=X.g;
    c=X.c_cell_type_tx;
    X=X.X;
end


[M,OUT]=run.i_talkr(X,g,c);

[n]=size(OUT.ligand_mat,2);
M=M.*log2(M*(n^2));
M(isnan(M))=0;
KL=real(sum(M,2));
Tok=[OUT.Tok, table(KL)];
[Tok,idx]=sortrows(Tok,'KL','descend');
%%
OUT.ligand_mat=OUT.ligand_mat(idx,:);
OUT.receptor_mat=OUT.receptor_mat(idx,:);
OUT.ligandok=OUT.ligandok(idx);
OUT.receptorok=OUT.receptorok(idx);
% assert(isequal(size(ligand_mat,1),length(ligandok)))
OUT.KL=KL(idx);
OUT.T=Tok;

end