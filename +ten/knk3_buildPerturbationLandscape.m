function [F]=knk3_buildPerturbationLandscape(A,genelist)
    import ten.*
n=length(genelist);
assert(n==size(A,1));
F=zeros(n,n);
    fprintf('\n');
    for k=1:n
        fprintf('Knocking out %s ... gene %d of %d\n',...
            genelist(k),k,n);
        [t]=ten.knk2_knockoutTargetGene(A,genelist,genelist(k),false);
        F(:,k)=t.drdist;
    end
end
