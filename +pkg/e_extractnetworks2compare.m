function [A0,A1,glist]=e_extractnetworks2compare(A0,A1,genelist,glist)
if nargin<4
    [glist]=gui.i_selectngenes(genelist);
end

glist=glist(:);
genelist=genelist(:);

[y,idx]=ismember(glist,genelist);
if any(~y), warning('Some genes in GLIST not found in GENELIST'); end
assert(isequal(genelist(idx(y)),glist(y)))

idx=idx(y);
A0=A0(idx,idx);
A1=A1(idx,idx);
glist=glist(y);

end
