function [Tup,Tdn]=e_processDETable(T,isatac)

if nargin<2, isatac=false; end

%a=T.pct_2-T.pct_1;
a=T.(T.Properties.VariableNames{8})-T.(T.Properties.VariableNames{7});
    if isatac
        isok=T.p_val_adj<=0.1;
    else
        isok=(abs(a)>=0.05 | abs(T.avg_log2FC)>=1.0) & T.p_val_adj<=0.01;
    end

    Tup=T(T.avg_log2FC > 0 & isok,:);
    Tdn=T(T.avg_log2FC < 0 & isok,:);
    
    Tup = sortrows(Tup,'abs_log2FC','descend');
    Tup = sortrows(Tup,'p_val_adj','ascend');

    Tdn = sortrows(Tdn,'abs_log2FC','descend');
    Tdn = sortrows(Tdn,'p_val_adj','ascend');
end