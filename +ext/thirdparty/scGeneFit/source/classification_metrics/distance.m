function [markers]= significant_cluster(data, labels, filename)
%data: [dxN]
%delta: [dxN]
close all
figure(1)

samples=double(data');
[n_total, d] = size(samples);
class_num = length(unique(labels));

% separate samples into different arrays based on label
separated = cell(1, class_num);
for i=1:n_total
    separated{labels(i)} = cat(1, separated{labels(i)}, samples(i,:));
end

markers=zeros(class_num,1);

bins=0:0.5:9;
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
    
    
    subplot(ceil(sqrt(class_num)),ceil(sqrt(class_num)),c)
    hold off
    histogram(current_class(:,markers(c)), bins, 'Normalization','probability')
    hold on;
    histogram(other_classes(:,markers(c)), bins, 'Normalization','probability')
    drawnow
end
saveas(gcf, filename)
end
