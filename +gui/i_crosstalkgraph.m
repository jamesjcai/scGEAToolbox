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
    

    figname=sprintf('%s (ligand) -> %s (receptor)',...
        ligandok(k),receptorok(k));    
    sc_grnview(m,cL,figname)
    
    figure;
    m=m-diag(diag(m));
    h=heatmap(m);
    
    h.XDisplayLabels=cL;
    h.YDisplayLabels=cL;    
    % pkg.heatmap(m, cL, cL,'%0.2f', 'TickAngle', 45, 'ShowAllTicks', true, 'TextColor', 'w');
end

