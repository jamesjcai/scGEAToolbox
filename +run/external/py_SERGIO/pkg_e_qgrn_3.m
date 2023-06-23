% C:\ProgramData\MATLAB\SupportPackages\R2023a\toolbox\matlab\quantum\

layer1=[];
for k=1:4
    layer1 = [layer1; ryGate(k,2*asin(sqrt(f0(k))))];
end

a=nchoosek(1:4,2);
theta0=rand(12,1);

layer2=[];
c=1;
for k=1:size(a,1)
    layer2=[layer2; cryGate(a(k,1),a(k,2),theta0(c))];
    c=c+1;
    layer2=[layer2; cryGate(a(k,2),a(k,1),theta0(c))];
    c=c+1;
end

C = quantumCircuit([layer1; layer2]);



np=12;
options = optimoptions('fmincon','Display','iter');
lb=-pi*ones(np,1);
ub=pi*ones(np,1);
x0=rand(np,1)*pi;

figure;
bar([ft, i_getp(C,x0)]);
[x,fval] = fmincon(@i_obj,x0,[],[],[],[],lb,ub,[],options,C,ft,layer1,layer2);

%options = optimset('Display','iter');
%[x,fval] = fminsearch(@i_obj,x0,options,C,ft,layer1,layer2);
%[x,fval] = fminunc(@i_obj,x0,options,C,ft);

figure;
bar([ft, i_getp(C,x)]);
[x0 x]'

i_getp(C,x)'
x(1:10)=pi./2;
i_getp(C,x)'

function [y]=i_obj(x,C,ft,layer1,layer2)
    C = quantumCircuit([layer1; layer2]);
    for k=(4+1):length(C.Gates)
        C.Gates(k).Angles=x(k-4);
    end
    S = simulate(C);
    [~,fo] = querystates(S);
    % f1=[probability(S,1,"1") probability(S,2,"1") ...
    %     probability(S,3,"1") probability(S,4,"1")];     % per gene activate freq.
    KL1 = sum(fo .* (log(fo)-log(ft)));
    KL2 = sum(ft .* (log(ft)-log(fo)));
    %y=mean([KL1 KL2])+norm(x);
    %[ft'; fo']
    y=100*mean([KL1 KL2]);
end
    

function [p]=i_getp(C,x)
    for k=(4+1):length(C.Gates)
        C.Gates(k).Angles=x(k-4);
    end
    S = simulate(C);
    [~,p] = querystates(S);
end
