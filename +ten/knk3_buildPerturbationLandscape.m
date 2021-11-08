function [F]=knk3_buildPerturbationLandscape(A,genelist,oldF,a,b)
if nargin<5, b=size(A,1); end
if nargin<4, a=1; end
if nargin<3, oldF=[]; end
import ten.*
n=length(genelist);
assert(n==size(A,1));
if ~isempty(oldF)
    F=oldF;
else
    F=zeros(n,n);
end
tmpmat=tempname;
    fprintf('Temporary result (tmpF) will be saved in %s.mat\n',tmpmat);
    a=min([a n]); b=min([b n]);
    for k=a:b
        fprintf('%s ... gene %d of %d\n',genelist(k),k,n);
        if all(F(:,k)==0)
            [t]=ten.knk2_knockoutTargetGene(A,genelist,genelist(k),false);
            F(:,k)=t.drdist;
        end
        if mod(k,100)==0
            try
                tmpF=F;
                fprintf('Saving temporary result (tmpF) to %s.mat\n',tmpmat);
                save(tmpmat,'tmpF','genelist');
            catch
            end
        end
    end
    try
        tmpF=F;
        fprintf('Saving temporary result (tmpF) to %s.mat\n',tmpmat);  
        save(tmpmat,'tmpF','genelist');
    catch
    end
end
