% C:\ProgramData\MATLAB\SupportPackages\R2023a\toolbox\matlab\quantum\


layer1=[];
for k=1:4, layer1 = [layer1; ryGate(k,2*asin(sqrt(f0(k))))]; end

%{
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
        [xa,fval] = fmincon(@i_obj,x0,[],[],[],[],lb,ub,[],options,pt,layer1);
    case 2
        options = optimset('Display','iter');
        [xa,fval] = fminsearch(@i_obj,x0,options,pt,layer1);
    case 3
        options = optimset('Display','iter');
        [xa,fval] = fminunc(@i_obj,x0,options,pt,layer1);
end


[poa]=i_fullcirc(xa,layer1);

subplot(2,2,4)
bar([pt poa])
set(gca,'XTick',1:length(states));
set(gca,'XTickLabel',states);
ylabel('# of cells');
xlabel('Expression pattern');
legend({'Target','Observed'})
title(sprintf('kl=%f',i_kldiverg(pt,poa)));


At=zeros(4);
c=1;
a=nchoosek(1:4,2);
for k=1:size(a,1)
    At(a(k,1),a(k,2))=xa(c);
    c=c+1;
    At(a(k,2),a(k,1))=xa(c);
    c=c+1;
end

figure;
subplot(2,2,1)
imagesc(abs(A)+abs(A'))
subplot(2,2,2)
imagesc(abs(At)+abs(At'))

% idx=eye(4);
% B=zeros(4);
% B(~idx)=xa;
% 
% c=1;
% xb=zeros(size(xa));
% 
% a=nchoosek(1:4,2);
% for k=1:size(a,1)
%     xb(c)=B(a(k,1),a(k,2));
%     c=c+1;
%     xb(c)=B(a(k,2),a(k,1));
%     c=c+1;
% end
% 


function [y]=i_obj(x,pt,layer1)
    [po,f1]=i_fullcirc(x,layer1);
    % f1=[probability(S,1,"1") probability(S,2,"1") ...
    %     probability(S,3,"1") probability(S,4,"1")];     % per gene activate freq.
    y1=i_kldiverg(pt,po,true);
    %[f0; f1]
    %pause
    f0=[0.4878    0.2954    0.4148    0.2834];
    y2=i_kldiverg(f0,f1,false);
    y=y1+y2;
end



