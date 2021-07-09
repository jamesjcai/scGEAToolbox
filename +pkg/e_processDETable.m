function [Tup,Tdn]=e_processDETable(T)

    isok=(abs(T.pct_2-T.pct_1)>0.05 | abs(T.avg_log2FC)>1.0)&T.p_val_adj<0.01;
    Tup=T(T.avg_log2FC>0 & isok,:);
    Tdn=T(T.avg_log2FC<0 & isok,:);
    
    Tup = sortrows(Tup,'abs_log2FC','descend');
    Tup = sortrows(Tup,'p_val_adj','ascend');

    Tdn = sortrows(Tdn,'abs_log2FC','descend');
    Tdn = sortrows(Tdn,'p_val_adj','ascend');
end