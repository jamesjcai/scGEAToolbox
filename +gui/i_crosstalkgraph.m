function i_crosstalkgraph(OUT,k)
if nargin<2, k=1; end
    cL=OUT.cL;
    ligand_mat=OUT.ligand_mat;
    receptor_mat=OUT.receptor_mat;
    ligandok=OUT.ligandok;
    receptorok=OUT.receptorok;
    
    a=ligand_mat(k,:);
    b=receptor_mat(k,:);
    m=(a'*b).*((a>0)'*(b>0));    
    m=m./sum(m(m>0));
    sc_grnview(m,cL)
    title(sprintf('%s (ligand) -> %s (receptor)',...
        ligandok(k),receptorok(k)));
end
