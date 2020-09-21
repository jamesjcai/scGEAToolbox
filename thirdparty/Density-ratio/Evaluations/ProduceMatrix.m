function [matrix,unassign,NetClass,noise] = ProduceMatrix (target,output)
%ProduceMatrix - generates a confusion matrix 
%
%IN:  target:   true labels
%     output:   clutering results
%OUT: matrix:   confusion matrix 
%     unassign:  number of unassigned instances
%     NetClass:  paried results for unassigned instances
%     noise:  paried results for unassigned instances

NetClass=[];
unassign=0;
matrix=zeros(max(target),max(output));
noise=zeros(max(target),1);
s=size(output,1);
for i = 1:s % test each example
    if output(i)>0
    matrix(target(i),output(i))=matrix(target(i),output(i))+1;
    NetClass=[NetClass;[target(i) output(i)]];
    else
        noise(target(i))=noise(target(i))+1;
    end
end
unassign=sum(noise);
% 
