function [ tpr ] = cal_specific_accuracy( W, labels, NUM )

[A,B] = sort(W + eye(length(W)),2,'descend');

if size(labels,1)<size(labels,2)
    labels = labels';
end

for i = 1 : NUM
ALLLABES(:,i) = labels(B(:,i));
trueindicator = ALLLABES(:,1:i) - repmat(labels,1,i);
tpr(i) = sum(trueindicator(:)==0)/(i*length(labels));
end

end

