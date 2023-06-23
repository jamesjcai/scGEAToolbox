function [Delta,smallest, lambda]=select_center_constraints(data, labels, samples_fraction, core)

[d,~]=size(data);
L=max(labels);
centers=zeros(d, L);

%compute the centers with all data
for i=1:L
    centers(:,i)=mean(data(:, labels==i), 2);
end

%randomly subsample data to produce constraints
[data, labels]=subsample(data, labels, samples_fraction );
[d,N]=size(data);

%constraints are of the form
%||Proj(point -other_center)|| - ||Proj(point -own_center)|| >= delta
smallest=intmax;
Delta=zeros(d, N*(L-1));
lambda=zeros(N*(L-1),1);
idx=1;
idx_prev=1;
for i=1:L
    aux=data(:, labels==i);
    for t=1:L
        if t~=i
            %iterate on all datapoints with label=i
            for j=1:size(aux,2)
                dist_own=aux(:,j)- centers(:,i);
                dist_other=aux(:,j)- centers(:,t);
                Delta(:,idx)=dist_other.*dist_other - 1.25* dist_own.*dist_own;
                
                s=sum(dist_own.*dist_own);
                if s<smallest
                    smallest=s;
                end
                idx=idx+1;
            end
            s=idx-idx_prev;
            lambda(idx_prev:idx-1)=s;
            idx_prev=idx;
        end
    end
end
end