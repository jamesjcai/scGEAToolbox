function [T]=e_grptest(y,grp)


[p_anova]=anova1(y,grp);
[p_kruskalwallis]=kruskalwallis(y,grp);
T=table(p_anova,p_kruskalwallis);
end