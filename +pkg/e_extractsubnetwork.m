function [outA,outglist]=e_extractsubnetwork(A,genelist,glist,stepk)

if nargin<4, stepk=0; end
glist=glist(:);
genelist=genelist(:);
[y,idxv]=ismember(glist,genelist);

outglist=glist(y);

if any(~y), warning('There are genes in GLIST not found in GENELIST'); end
assert(isequal(genelist(idxv(y)),outglist))

idxv=idxv(y);
if stepk>0
    for k=1:length(idxv)
        a=abs(A(idxv(k),:));
        b=abs(A(:,idxv(k)));
        [~,mid]=maxk(a,stepk);
        outglist=[outglist; genelist(mid)];
        [~,mid]=maxk(b,stepk);
        outglist=[outglist; genelist(mid)];
        outglist=unique(outglist);
    end
    [~,idxv]=ismember(outglist,genelist);
end
    outA=A(idxv,idxv);
end
