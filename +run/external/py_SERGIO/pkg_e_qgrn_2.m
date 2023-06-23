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
% C = quantumCircuit([layer1]);
%figure; plot(C);


S = simulate(C);
%S.BasisStates
%S.Amplitudes
% f = formula(S)
% figure; histogram(S)

[states,P] = querystates(S);
assert(isequal(states,S.BasisStates))
assert(isequal((S.Amplitudes).^2 ,P))
assert(isequal(txt,states))

figure; bar(P)
set(gca,'XTick',1:length(states));
set(gca,'XTickLabel',states);
ylabel('# of cells');
xlabel('Expression pattern');

f1=[probability(S,1,"1") probability(S,2,"1") ...
    probability(S,3,"1") probability(S,4,"1")];     % per gene activate freq.
% figure; bar([f0; f1]')
fo=P;


%M = randsample(S,50)
%T = table(M.Counts,M.MeasuredStates,VariableNames=["Counts","States"])




