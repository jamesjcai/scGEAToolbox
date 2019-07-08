% function [NMI,ARI,grps,similarity] = SinNLRR(path_data,path_label)
% 
% in_X=load(path_data);
% label=load(path_label);
% n_space = length(unique(label));
% [X,] = FilterGenesZero(in_X);
% [n,m]=size(X);
% X = normalize(X');
function [grps] = SinNLRR(X,k)


% [n,m]=size(X);
% X = normalize(X');
% 
% in_X=load('single_data/Test_12_Goolam.txt');
% label=load('single_data/Test_12_Goolam_label.txt');
% [X,] = FilterGenesZero(in_X);
% [n,m]=size(X);
% X = normalize(X');

X=X./vecnorm(X);
[n]=size(X,2);
n_space=k;


% NMI=0;
% ARI=0;
    for i=1:15
        alpha=0.1+i*0.2;
        Z = lrr_relaxed(X,alpha);
        [localX, ~] = localize(abs(Z));
        %z_rate = sum(sum(localX~=0))/(n*n);
        min_neighbour = min(sum(localX~=0,1));
        %max_neighbour = max(sum(localX~=0,1));
        if ((n<1000 &&min_neighbour>3)||(n>=1000 &&min_neighbour>1))          
            similarity = localX+localX';    
            grps = SpectralClustering(similarity, n_space);
%             NMI=Cal_NMI(label, grps);
%             ARI=Cal_ARI(label, grps);
            break;
        end
        %fprintf('%f\t%f\t%f\t%f\t%f\t%f\t%f\n',alpha,z_rate,min_neighbour,max_neighbour,mean(col_cov),missrate,misARI);
    end
end