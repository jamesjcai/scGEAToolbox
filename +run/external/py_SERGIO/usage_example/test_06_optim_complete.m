g2=ryGate(1:size(a2,1),a2(:,1));
c2=quantumCircuit(g2);
s2=simulate(c2);
[states2,p2] = querystates(s2,1:10,Threshold=0);


% fun = @(x)3*x(1)^2 + 2*x(1)*x(2) + x(2)^2 - 4*x(1) + 5*x(2);
% x0 = [1,1];
% [x,fval] = fminunc(fun,x0)
x0=rand(90,1)*pi;
x0=zeros(90,1);
g1=ryGate(1:size(a1,1),a1(:,1));
options = optimset('Display','iter');  %'MaxIter', 200*90

%[x,fval] = fminunc(@i_obj,x0,options,g1,p2);

options = optimoptions('fmincon','Display','iter');
lb=-pi*ones(90,1);
ub=pi*ones(90,1);
[x,fval] = fmincon(@i_obj,x0,[],[],[],[],lb,ub,[],options,g1,p2);


% n1=6; n2=1;
% (n1-1)*9+n2

%[4 6]
%[4 2]
%[5 3]
ff = @(x)(x(1)-1)*9+x(2)-(x(2)>x(1));
idx=[ff([4 6]) ff([6 4]) ff([4 2]) ff([2 4]) ff([5 3]) ff([3 5])];
[x(idx) x(idx+1)]
figure;
plot(x);
xline(idx,'r');

function y=i_obj(x,g1,p2)
    
    gates=g1;
    c=0;
    for k1=1:10
        for k2=1:10
            if k1~=k2
                c=c+1;
                gx=cryGate(k1,k2,x(c));
                gates=[gates;gx];
            end
        end
    end

    c1=quantumCircuit(gates);
    s1=simulate(c1);
    [~,p1] = querystates(s1,1:10,Threshold=0);
    KL1 = sum(p1 .* (log(p1)-log(p2)));
    KL2 = sum(p2 .* (log(p2)-log(p1)));
    y=mean([KL1 KL2])+0*norm(x);
end


