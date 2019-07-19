function [ total_data_out] = SC_reduce( total_data, blocks )
% SC_reduce ************************************************************************
% GIOVANNI IACONO, CNAG, 16/08/2017
% generates convoluted expression matrix give an an expression matrix and
% the indexes of its iCells

% INPUT
% total_data: expression matrix, not normalized
% blocks: indices of iCella

% OUTPUT
% total_data_out: convoluted expression matrix


num_genes=length(total_data(:,1));
num_samples=length(total_data(1,:));

total_data_out=zeros(num_genes,length(blocks));

for k=1:length(blocks)
total_data_out(:,k)=sum(total_data(:,nonzeros(blocks(k,:)))');
end




end