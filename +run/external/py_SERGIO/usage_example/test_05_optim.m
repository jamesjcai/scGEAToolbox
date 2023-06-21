g2=ryGate(1:size(a2,1),a2(:,1));
c2=quantumCircuit(g2);
s2=simulate(c2);
[states2,p2] = querystates(s2,1:10,Threshold=0);


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

% fun = @(x)3*x(1)^2 + 2*x(1)*x(2) + x(2)^2 - 4*x(1) + 5*x(2);
% x0 = [1,1];
% [x,fval] = fminunc(fun,x0)
x0=rand(90,1)*pi;
x0=zeros(90,1);
options = optimset('Display','iter');
[x,fval] = fminsearch(@i_obj,x0,options,gates,p2);

% n1=6; n2=1;
% (n1-1)*9+n2

function y=i_obj(x,gates,p2)
    
    %y=3*x(1)^2 + 2*x(1)*x(2) + x(2)^2 - 4*x(1) + 5*x(2);
    for k=1:90
        gates(k+10,:).Angles=x(k);
    end
    gates(10:14,:)
    c1=quantumCircuit(gates);
    s1=simulate(c1);
    [~,p1] = querystates(s1,1:10,Threshold=0);
    KL1 = sum(p1 .* (log(p1)-log(p2)));
    KL2 = sum(p2 .* (log(p2)-log(p1)));
    y=mean([KL1 KL2]);
end


