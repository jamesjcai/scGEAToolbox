function [po,f1]=i_fullcirc(theta0,layer1)
    n=size(layer1,1);
%    layer1=[];
%    for k=1:4, layer1 = [layer1; ryGate(k,2*asin(sqrt(f0(k))))]; end
    a=nchoosek(1:n,2);
    layer2=[];
    c=1;
    for k=1:size(a,1)
        layer2=[layer2; cryGate(a(k,1),a(k,2),theta0(c))];
        c=c+1;
        layer2=[layer2; cryGate(a(k,2),a(k,1),theta0(c))];
        c=c+1;
    end
    C = quantumCircuit([layer1; layer2]);
    S = simulate(C);
    [~,po] = querystates(S);
    if nargout>1        
        f1=zeros(1,n);
        for k=1:n
            f1(k)=probability(S,k,"1");
        end
    end
end