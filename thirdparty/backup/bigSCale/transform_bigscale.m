function [X]=transform_bigscale(X)
% model=1. Log(x), then each row (gene) normalized between [-5:5]
%
% adapted from:
%
% SC_log_transform ************************************************************************
% GIOVANNI IACONO, CNAG, 16/08/2017
% Transfomes a normalized expression matrix into different gene-normalized
% matrices
% INPUT
% indices: optional, restrict to "indices" genes. If empty uses all genes (Default).
% total_data_norm: normalized expression matrix
% model: see beloe for the descriptions
% OUTPUT 
% matrix: list of all significant markers
    [X]=norm_libsize(X);
    X=log2(X+1);
    X=10*(X-min(X,[],2))./(max(X,[],2)-min(X,[],2))-5;
    X(isnan(X))=0;
    
%     return;
%     for k=1:size(X,1)
%         if length(unique(X(k,:)))>2
%             x_out=shift_values(X(k,:),0,10);
%         else
%             x_out=zeros(size(X(k,:)));
%         end
%         X(k,:)=x_out;
%     end
% 	X=X-5;
% end
% 
% function x_out=shift_values(x,X,Y)
% a=min(x);
% b=max(x);
% x_out=((x-a)/(b-a))*(Y-X)+X;
% end 

