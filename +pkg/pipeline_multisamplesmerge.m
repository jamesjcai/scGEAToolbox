% https://www.cell.com/cell/fulltext/S0092-8674(19)31178-X
% DOI:https://doi.org/10.1016/j.cell.2019.10.028

%{'T11-Apobec-7daytreated','T11-Apobec-Nottreated',...
% 'KPB25Luv-7daytreated','KPB25Luv-Nottreated'}

accv={'GSM4042585','GSM4042586','GSM4042587','GSM4042588'};
SCEV=cell(length(accv),1);
for k=4:-1:3 % :length(accv)
    pause(5);
    [sce]=sc_readgeoaccession(accv{k});
    sce=sce.qcfilterwhitelist(500,0.15,0.01,'');
    SCEV{k}=sce;
end
%%
sce=sc_mergesces(SCEV);
sce = sce.embedcells('tsne', true, true, 3);
k=round(sce.NumCells/100);
sce = sce.clustercells(k, 'kmeans', true);
[sce]=pkg.e_celltypes2allclust(sce,'mouse',true);


