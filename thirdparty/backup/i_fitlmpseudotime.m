%load S2_data_G1_only.mat

n=length(genelist);
% n=5;
beta1=zeros(n,1);
p1=zeros(n,1);
corre=zeros(n,1);
p2=zeros(n,1);
for k=1:n
    k
    md=fitlm(t_mono,Ximp(k,:)');
    [r,p]=corr(t_mono,Ximp(k,:)');
    corre(k)=r;
    p2(k)=p;
    beta1(k)=md.Coefficients.Estimate(2);
    p1(k)=md.Coefficients.pValue(2);    
end
fdr1 = mafdr(p1);
fdr2 = mafdr(p2);
T=table(genelist,beta1,p1,fdr1,corre,p2,fdr2);
T=sortrows(T,2,'descend');

writetable(T,'pseudotime_res.txt');
