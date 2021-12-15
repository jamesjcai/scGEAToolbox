accv={'GSM4042585','GSM4042586','GSM4042587','GSM4042588'};
SCEV=cell(length(accv),1);
for k=1:length(accv)
    pause(5);
    [sce]=sc_readgeoaccession(accv{k});
    sce=sce.qcfilterwhitelist(500,0.15,0.01,'');
    SCEV{k}=sce;
end
sce=sc_mergesces(SCEV);
sce = sce.embedcells('tsne', true, true, 3);
k=round(sce.NumCells/100);
sce = sce.clustercells(k, 'kmeans', true);
[sce]=pkg.e_celltypes2allclust(sce,'mouse',true);


