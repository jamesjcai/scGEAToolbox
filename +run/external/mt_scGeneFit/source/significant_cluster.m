%Filename:      significant_cluster.m
%Last edit:     Feb 11 2019
%Description:   Implementation of the state of the art method from [1].
%Inputs:
%               -data: 
%     
%               a d x N array (d dimension, N number of points)
% 
%               -labels: 
% 
%               a 1 x N array with labels from 1 to L
%
%
% 
%Outputs:
%               -markers: 
% 
%               a 1 x L array with the indices of the selected markers
%
% 
%Documentation: [1]

function [markers]= significant_cluster(data, labels)


samples=double(data');
[n_total, d] = size(samples);
class_num = length(unique(labels));

% separate samples into different arrays based on label
separated = cell(1, class_num);
for i=1:n_total
    separated{labels(i)} = cat(1, separated{labels(i)}, samples(i,:));
end

markers=zeros(class_num,1);

bins=0:max(max(data))/20:max(max(data));
for c=1:class_num
    current_class = separated{c};
    other_classes = cat(1, separated{1:end ~= c});
    bigdist=0;
    for marker=1:d    
        [c1,n1]=hist(current_class(:,marker),bins);
        c1=c1/size(current_class(:,marker),1);
        [c2,n2]=hist(other_classes(:,marker),bins);
        c2=c2/size(other_classes(:,marker),1);
        
        dist=chi_square_statistics(c1,c2);
        if dist>bigdist
            bigdist=dist;
            markers(c)=marker;
        end
    end    
end

end
