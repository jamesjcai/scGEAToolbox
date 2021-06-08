function sc_qcviolin(X,genelist)

i=startsWith(genelist,'mt-','IgnoreCase',true);
nftr=sum(X>0,1);
lbsz=sum(X,1);
lbsz_mt=sum(X(i,:),1);
cj=100*(lbsz_mt./lbsz);
figure;
subplot(1,3,1)
pkg.violinplot(nftr,[],'showdata',false);
title('nFeatur\_RNA')
subplot(1,3,2)
pkg.violinplot(lbsz,[],'showdata',false);
title('nCount\_RNA')
subplot(1,3,3)
pkg.violinplot(cj,[],'showdata',false);
title('percent.mt')



