function [cL]=e_sortcatbysize(thisc)
    [c, cL] = grp2idx(thisc);
    cmv = 1:max(c);
    %idxx = cmv;
    [cmx] = countmember(cmv, c);
    [~, idxx] = sort(cmx, 'descend');
    cL=cL(idxx);
end
