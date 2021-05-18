function [outA0,outA1,outglist]=e_extractnetworks2compare(A0,A1,genelist,glist,stepk)
if nargin<5, stepk=0; end
if nargin<4
    [glist]=gui.i_selectngenes(genelist);
end

glist=glist(:);
genelist=genelist(:);

[y,idxv]=ismember(glist,genelist);

outglist=glist(y);

if any(~y), warning('There are genes in GLIST not found in GENELIST'); end
assert(isequal(genelist(idxv(y)),outglist))



% idx=idx(y);
% A0=A0(idx,idx);
% A1=A1(idx,idx);
% glist=glist(y);

idxv=idxv(y);
if stepk>0
    for k=1:length(idxv)
        a=abs(A0(idxv(k),:));
        [~,mid]=maxk(a,stepk);
        outglist=[outglist; genelist(mid)];
        b=abs(A0(:,idxv(k)));
        [~,mid]=maxk(b,stepk);
        outglist=[outglist; genelist(mid)];
        outglist=unique(outglist);

        a=abs(A1(idxv(k),:));
        [~,mid]=maxk(a,stepk);
        outglist=[outglist; genelist(mid)];
        b=abs(A1(:,idxv(k)));
        [~,mid]=maxk(b,stepk);
        outglist=[outglist; genelist(mid)];
        outglist=unique(outglist);        
    end
    [~,idxv]=ismember(outglist,genelist);
end
outA0=A0(idxv,idxv);
outA1=A1(idxv,idxv);
end
