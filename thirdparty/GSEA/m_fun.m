function ranks = m_fun(data,group)
% External function for calculating gene ranking metric. 
% Input:
% data - data matrix
% group - grouping sample
% Output:
% ranks - ranking metric values
% Author:
% Michal.Marczyk@polsl.pl

%Here you can put your own code for calculating gene ranking metric.
%For example:
ranks = mean(data(:,group==0),2) - mean(data(:,group==1),2);
