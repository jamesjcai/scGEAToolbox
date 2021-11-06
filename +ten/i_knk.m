function [T]=i_knk(A0,idx,genelist)

    if ischar(idx)||isstring(idx)
        [~,idx]=ismember(idx,genelist);
    end
    
    if sum(A0(idx,:)~=0)<50
        warning('KO gene (%s) has no link or too few links with other genes.',...
                 genelist(idx));
        T=table();
        return;
    end
    import ten.*
    A1=A0;
    A1(idx,:)=0;
    
    a=A0(idx,:); b=A1(idx,:);
    yn = abs(a - b) <= eps(max(abs(a), abs(b)));
    yn = all(yn);  % scalar logical output
    
    if yn
        drdist=nan(length(genelist),1);
        T=table(drdist);
    else
        [aln0,aln1]=i_ma(A0,A1);
        T=i_dr(aln0,aln1,genelist,true);
    end
end