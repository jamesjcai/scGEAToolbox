function [ci,cLi]=i_grp2idxsorted(s)

    [ci, cLi] = grp2idx(s);
    %---------
    [cLisorted,idx]=natsort(cLi);
    cisorted=ci;
    for k=1:length(idx), cisorted(ci==idx(k))=k; end
    ci=cisorted;
    cLi=cLisorted;

end
