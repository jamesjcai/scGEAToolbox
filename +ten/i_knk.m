function [T]=i_knk(A0,idx,genelist,dosort,lambdav)

%see also: ten.knk2_knockoutTargetGene

    if nargin<5, lambdav=0; end
    if nargin<4, dosort=true; end
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

    if lambdav~=0
        S=A0.*(abs(A0)>abs(A0.'));
        A0=(1-lambdav)*A0+lambdav*S;
    end
    A0=A0-diag(diag(A0));
    A0=A0.';
    
    A1=A0;
    A1(idx,:)=0;
    
    a=A0(idx,:); b=A1(idx,:);
    yn = abs(a - b) <= eps(max(abs(a), abs(b)));
    yn = all(yn);     % scalar logical output
    
    if yn
        drdist=nan(length(genelist),1);
        T=table(drdist);
    else
        [aln0,aln1]=i_ma(A0,A1);
        T=ten.i_dr(aln0,aln1,genelist,dosort);
    end
end