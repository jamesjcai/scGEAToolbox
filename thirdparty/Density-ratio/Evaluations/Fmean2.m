function [ fmeasure,recall,precision ] = Fmean2( matrix,match )
% FMEAN Computes the fmeasure,recall,precision based on the confusion matrix (matrix)
% match is the best matching on the confusion matrix

%recall
[r,c]=size(matrix);
[rr,cc]=size(match);
match=[match;zeros(r-rr,cc)];
match=[match zeros(r,c-cc)];


recall=matrix;
sumrow=sum(matrix')';
if size(matrix,2)==1
    sumrow=matrix;
end
for j = 1:size(matrix,1)
        recall(j,:) = recall(j,:)/(sumrow(j)+0.0000001); % calculate the total positive examples 
end

%precision
precision=matrix;
sumcol=sum(matrix);
if size(matrix,1)==1
    sumcol=matrix;
end
for j = 1:size(matrix,2)
        precision(:,j) = precision(:,j)/(sumcol(j)+0.0000001); % calculate the total positive examples 
end
%fmeasure

fmeasure=2*precision.*recall./(precision+recall+0.0000001);
fmeasure=fmeasure(match==1);
precision=precision(match==1);
recall=recall(match==1);
end

