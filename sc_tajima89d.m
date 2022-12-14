%targetg=
%g


idx=[1 2 3 4 5 6 7 8];
x=sce.X(idx,:);
g=sce.g(idx);

[smpln]=size(x,2);
p=1-sum(x>0,2)./smpln;
Sn=sum(p>0 & p<1);

x1=sum(2.*p.*(1-p));
x2=sum(2.*p.*p);
y=smpln/(smpln-1);
thepi=x1*y;
theh=x2*y;

[h]=thepi-theh;
% [d,thetw,pval] = tajima89d(smpln, Sn, thepi)
% [d] = tajima89d(smpln, Sn, thepi)
nx=1:(smpln-1);
a1 = sum(1./nx);
[d]=(thepi-Sn/a1);


