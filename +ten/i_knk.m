function [T]=i_knk(A0,idx,genelist)

    if ischar(idx)||isstring(idx)
        [~,idx]=ismember(idx,genelist);
    end

    import ten.*
    A1=A0;
    A1(idx,:)=0;
    [aln0,aln1]=i_ma(A0,A1);
    T=i_dr(aln0,aln1,genelist,true);

end