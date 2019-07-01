function [x,sf]=norm_deseq(x)
% For DESeq normalization, the geometric mean for each gene was computed 
% after removing all zeroes. This is necessary to avoid a situation where 
% a majority of genes have geometric means of zero, such that the majority
% of ratios to the geometric mean would be undefined. Size factors were 
% then computed using the estimateSizeFactorsForMatrix function in DESeq2 
% v1.10.1 [21]. In this function, ratios of zero were automatically removed
% prior to calculation of the median in each library, to avoid obtaining a
% size factor equal to zero        
    x(x==0)=nan;
    sf=nanmedian(x./nangeomean(x,2));
    x=x./sf;
    x(isnan(x))=0;
end


