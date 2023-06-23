% C:\ProgramData\MATLAB\SupportPackages\R2023a\toolbox\matlab\quantum\

%{
layer1=[];
for k=1:4, layer1 = [layer1; ryGate(k,2*asin(sqrt(f0(k))))]; end


theta0=i_randpmpi(12,1);
a=nchoosek(1:4,2);
layer2=[];
c=1;
for k=1:size(a,1)
    layer2=[layer2; cryGate(a(k,1),a(k,2),theta0(c))];
    c=c+1;
    layer2=[layer2; cryGate(a(k,2),a(k,1),theta0(c))];
    c=c+1;
end
C = quantumCircuit([layer1; layer2]);

theta0=i_randpmpi(12,1)';
[po]=i_qrstates(C,theta0)';
theta0=i_randpmpi(12,1)';
[po]=i_qrstates(C,theta0)';
theta0=i_randpmpi(12,1)';
[po]=i_qrstates(C,theta0)';

return;

%}

%[pt po]'
%[kl]=i_kldiverg(pt,po)


methodid=3;
np=12;
x0=i_randpmpi(np,1);

switch methodid
    case 1        
        options = optimoptions('fmincon','Display','iter');
        lb=-pi*ones(np,1);
        ub=pi*ones(np,1);
        [xa,fval] = fmincon(@i_obj,x0,[],[],[],[],lb,ub,[],options,pt,f0);
    case 2
        options = optimset('Display','iter');
        [xa,fval] = fminsearch(@i_obj,x0,options,pt,f0);
    case 3
        options = optimset('Display','iter')
        [xa,fval] = fminunc(@i_obj,x0,options,pt,f0);
end


[poa]=i_fullcirc(xa,f0);

subplot(2,2,4)
bar([pt poa])
set(gca,'XTick',1:length(states));
set(gca,'XTickLabel',states);
ylabel('# of cells');
xlabel('Expression pattern');
legend({'Target','Observed'})
title(sprintf('kl=%f',i_kldiverg(pt,poa)));

function [y]=i_obj(x,pt,f0)
    [po]=i_fullcirc(x,f0);
    % f1=[probability(S,1,"1") probability(S,2,"1") ...
    %     probability(S,3,"1") probability(S,4,"1")];     % per gene activate freq.
    y=i_kldiverg(pt,po);
end



