function [out_X,flag_index] = FilterGenesZero(in_X)
    [m,n]=size(in_X);
    flag=(in_X~=0);
    flag_count=sum(flag,1);
    flag_index=flag_count>(m*0.05);
    out_X=in_X(:,flag_index);
    %out_Genes=Genes(:,flag_count>300);
end