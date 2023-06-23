% C:\ProgramData\MATLAB\SupportPackages\R2023a\toolbox\matlab\quantum\


% f0 is used to initialize the first layer ry gates.
layer1=[];
for k=1:4, layer1 = [layer1; ryGate(k,2*asin(sqrt(f0(k))))]; end


%theta0=rand(12,1);
theta0=pi*((rand(12,1)*2)-1);

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
% C = quantumCircuit([layer1]);
%figure; plot(C);


S = simulate(C);
%S.BasisStates
%S.Amplitudes
% f = formula(S)
% figure; histogram(S)

[states,po] = querystates(S);
assert(isequal(states,S.BasisStates))
assert(isequal((S.Amplitudes).^2, po))
assert(isequal(txt(:),states(:)))

%figure; 
subplot(2,2,2)
bar(po)
set(gca,'XTick',1:length(states));
set(gca,'XTickLabel',states);
ylabel('# of cells');
xlabel('Expression pattern');

%f1=[probability(S,1,"1") probability(S,2,"1") ...
%    probability(S,3,"1") probability(S,4,"1")];     % per gene activate freq.
% figure; bar([f0; f1]')

%figure; 
subplot(2,2,3)
bar([pt po])
set(gca,'XTick',1:length(states));
set(gca,'XTickLabel',states);
ylabel('# of cells');
xlabel('Expression pattern');
legend({'Target','Observed'})
title(sprintf('kl=%f',i_kldiverg(pt,po)));


%M = randsample(S,50)
%T = table(M.Counts,M.MeasuredStates,VariableNames=["Counts","States"])




