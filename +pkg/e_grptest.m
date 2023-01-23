function [T]=e_grptest(y,grp)

    T=[];
    if length(unique(grp))==2
        p_ttest=zeros(size(y,1),1);
        p_wilcoxon=zeros(size(y,1),1);  % Mann Whitney U Test (Wilcoxon Rank Sum Test)
        id=grp2idx(grp);
        for k=1:size(y,1)
            a=y(k,id==1); b=y(k,id==2);
            [~,p]=ttest2(a,b);
            p_ttest(k)=p;
            p_wilcoxon(k)=ranksum(a,b);
            % The Kruskal-Wallis test is a nonparametric version of classical one-way ANOVA, and an extension of the Wilcoxon rank sum test to more than two groups.            
        end
        T=table(p_ttest,p_wilcoxon);    
    elseif length(unique(grp))>2
        p_anova=zeros(size(y,1),1);
        p_kruskalwallis=zeros(size(y,1),1);
        for k=1:size(y,1)
            p_anova(k)=anova1(y(k,:),grp,'off');
            p_kruskalwallis(k)=kruskalwallis(y(k,:),grp,'off');
            % The Kruskal-Wallis test is a nonparametric version of classical one-way ANOVA, and an extension of the Wilcoxon rank sum test to more than two groups.            
        end
        T=table(p_anova,p_kruskalwallis);
    end

end