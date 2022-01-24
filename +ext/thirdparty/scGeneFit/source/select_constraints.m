%Filename:      select_constraints.m
%Last edit:     Feb 11 2019
%Description:   For each point in the dataset we consider the k-nearest 
%               from each of the other clusters and we return the difference
%               vector
%Inputs:
%               -data: 
%     
%               a d x N array (d dimension, N number of points)
% 
%               -labels: 
% 
%               a 1 x N array with labels from 1 to L
%
%               -k: 
%     
%               the number of neighbors to select. If k=0 all constraints
%               are returned
%
%               -opt: 
%     
%               struct of optional parameters:
%               opt.samples_fraction (default 1): 
%                   fraction of data to get the constraints from.
%                   if smaller than 1 the data is selected at random.
%               opt.constraints_neighbors (default 4)
%                   number of nearest neighbors to select constraints from
%               opt.hinge_scale (default 0.1)
%                   hinge parameter for linear program (see paper)
% 
%Outputs:
%               -Delta: 
% 
%               a d x n array, where each column is of the form (v-w).^2  
%               where v and w are nearest neighbors with different labels
%
%               -smallest: 
% 
%               the smallest norm of vectors in Delta 
% 
%Documentation:

function [Delta, smallest]= select_constraints(data, labels, k)

samples=double(data');
[n_total, d] = size(samples);
class_num = length(unique(labels));

% separate samples into different arrays based on label
separated = cell(1, class_num);
for i=1:n_total
    separated{labels(i)} = cat(1, separated{labels(i)}, samples(i,:));
end


idx=1;
smallest=intmax;

if k==0
    Delta=zeros(d,0);
    %all constraints
    for i=1:n_total
        for j=i+1:n_total
            if labels(i)~=labels(j)
                delta=data(:,i)-data(:,j);
                Delta(:,idx)=delta;
                if smallest>norm(delta)^2
                    smallest=norm(delta)^2;
                end
                idx=idx+1;
            end
        end
    end
else
    Delta=zeros(d,k*n_total);
    
    for c=1:class_num
        % Separate data by current and other classes
        current_class = separated{c};
        for cc=c+1:class_num
            %maybe change this to cc=1:class_num
            if cc~=c
                other_classes = separated{cc};
                
                % Search for kNN
                nn_idxs = knnsearch(other_classes, current_class, 'k', k);
                
                for i=1:size(current_class,1)
                    % Loop over number of nearest neighbors
                    for j=1:k
                        % Get current nearest neighbor index in nn_idxs
                        nn_idx = nn_idxs(i, j);
                        delta = current_class(i,:) - other_classes(nn_idx, :);
                        if smallest>norm(delta)^2
                            smallest=norm(delta)^2;
                        end
                        Delta(:,idx)=double(delta');
                        idx=idx+1;
                    end
                end
            end
        end
    end
    Delta=Delta.*Delta;
end
