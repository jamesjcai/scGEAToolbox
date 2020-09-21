function [ matrix,unassign,match,nmiscore,accurateScore,NetFmeasure,recall,precision,Fmeasure,indF] = evaluate(class,Tclass)
%EVALUATE computes different evaluation measure based on the true label and cluster label
%   Input    : true_labels    : N-by-1 vector containing true labels
%              cluster_labels : N-by-1 vector containing cluster labels
%   Output   : matrix: confusion matrix
%              unassign : the number of unassigned instances
%              match : best match between the true label and cluster label
%              nmiscore: normalised mutual information (NMI) 
%              accurateScore : accurate score
%              NetFmeasure : F1 measure without unassigned instances
%              recall : recall score
%              precision : precision score
%              Fmeasure : F1 measure with unassigned instances
%              indF: individual F1 measure (with unassigned instances) for each cluster

[matrix,unassign,NetClass,noise] = ProduceMatrix (class,Tclass);
indF=zeros(max(class),1);
if size(matrix,2)~=0
    nmiscore=nmi(NetClass(:,1), NetClass(:,2));
    [accurateScore, ~] = accuracy(NetClass(:,1), NetClass(:,2));
    accurateScore=accurateScore/100;
    
    [nFmeasure,~,~,match]=Fmean(matrix);
    NetFmeasure=sum(nFmeasure)/max(class);
    
    
    [Fmeasure,recall,precision]=Fmean2( [matrix noise],match);
    
    %
    [posc,posr]=find(match==1);
    for i=1:size(posc,1)
        indF(posc(i))=Fmeasure(i);
    end
    
    Fmeasure=sum(indF)/max(class);
    recall=sum(recall)/max(class);
    precision=sum(precision)/max(class);
    
else
    match=0;
    nmiscore=0;
    accurateScore=0;
    NetFmeasure=0;
    Fmeasure=0;
    recall=0;
    precision=0;
end

end

