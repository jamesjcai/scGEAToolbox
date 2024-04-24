function [A0, A1, glist] = i_extractnet2compare(A0, A1, genelist, glist)
if nargin < 4
    [glist] = gui.i_selectngenes(genelist);
end
[y, idx] = ismember(glist, genelist);
if any(~y), error('xx'); end
A0 = A0(idx, idx);
A1 = A1(idx, idx);
end