function gui_qcfilter(X,genelist)

figure; 
nGenes=sum(X>0);
nUMIs=sum(X);
mtratio=sc_mtratio(X,genelist);

subplot(2,2,1)
scatter(nGenes,mtratio);
xlabel('nGenes'); ylabel('pMito');
xline(300,'r--'); xline(5000,'r--');
yline(0.025,'r--'); yline(0.1,'r--');

subplot(2,2,2)
scatter(nUMIs,mtratio);
xlabel('nUMI'); ylabel('pMito');
xline(400,'r--'); xline(30000,'r--');
yline(0.025,'r--'); yline(0.1,'r--');


subplot(2,2,3)
scatter(nGenes,nUMIs);
xlabel('nGenes'); ylabel('nUMI');
xline(300,'r--'); xline(5000,'r--');
yline(400,'r--'); yline(30000,'r--');

% https://github.com/prabhakarlab/RCAv2