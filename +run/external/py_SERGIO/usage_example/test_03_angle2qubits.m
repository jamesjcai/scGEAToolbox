%tan(c1(:,2)./c1(:,1))
%cart2sph(c1(:,1),c1(:,2),0)

g1=ryGate(1:size(a1,1),a1(:,1));
c1=quantumCircuit(g1);
%plot(c)
s1=simulate(c1);
[states1,p1] = querystates(s1,1:10,Threshold=0);
%close all; 

g2=ryGate(1:size(a2,1),a2(:,1));
c2=quantumCircuit(g2);
%plot(c)
s2=simulate(c2);
[states2,p2] = querystates(s2,1:10,Threshold=0);
%close all; 

% figure;
% subplot(1,2,1)
% histogram(s1)
% subplot(1,2,2)
% histogram(s2)

%%
KL1 = sum(p1 .* (log(p1)-log(p2)));
KL2 = sum(p2 .* (log(p2)-log(p1)));

(KL1+KL2)/2

%%
% KL_div = sum(p1 .* log2(p1 ./ p2), 'omitnan')
% 
% P = [.05, .1, .2, .05, .15, .25, .08, .12]
% Q = [.3, .1, .2, .1, .1, .02, .08, .1]
% KL = sum(P .* (log2(P)-log2(Q)));
% KL = sum(P .* (log(P)-log(Q)))
% 
