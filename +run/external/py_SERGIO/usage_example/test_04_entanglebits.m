g1=ryGate(1:size(a1,1),a1(:,1));
gates=g1;
for k1=1:10
    for k2=1:10
        if k1~=k2
            gx=cryGate(k1,k2,0);
            gates=[gates;gx];
        end
    end
end

     for k=1:90
         gates(k+10).Angles=rand;
     end
    % gates(100).Angles=3;

%n=10;
%n*(n-1)

%%
c1=quantumCircuit(gates);
figure;
plot(c1)

%%
s1=simulate(c1);
[states1,p1] = querystates(s1,1:10,Threshold=0);


%%

g2=ryGate(1:size(a2,1),a2(:,1));
c2=quantumCircuit(g2);
%plot(c)
s2=simulate(c2);
[states2,p2] = querystates(s2,1:10,Threshold=0);


KL1 = sum(p1 .* (log(p1)-log(p2)));
KL2 = sum(p2 .* (log(p2)-log(p1)));

KL=mean([KL1 KL2])

