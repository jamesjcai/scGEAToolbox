function i_crosstalkgraph2(OUT1,OUT2,k)
if nargin<3, k=1; end
    m1=i_getm(OUT1,k);
    m2=i_getm(OUT2,k);
    cL1=OUT1.cL;
    cL2=OUT2.cL;
    assert(isequal(cL1,cL2))
    sc_grnview2(m1,m2,cL1)
    title(sprintf('%s (ligand) -> %s (receptor)',...
        OUT1.ligandok(k),OUT1.receptorok(k)));
end


function m=i_getm(OUT,k)
    %cL1=OUT.cL;
    ligand_mat=OUT.ligand_mat;
    receptor_mat=OUT.receptor_mat;
    %ligandok=OUT.ligandok;
    %receptorok=OUT.receptorok;
    
    a=ligand_mat(k,:);
    b=receptor_mat(k,:);
    m=(a'*b).*((a>0)'*(b>0));
    m=m./sum(m(m>0));
end
