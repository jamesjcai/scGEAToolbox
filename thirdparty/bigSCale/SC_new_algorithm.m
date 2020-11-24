function [N, N_pct, matrix] = SC_new_algorithm( total_data , edges, ds)
% SC_new_algorithm ************************************************************************
% GIOVANNI IACONO, CNAG, 16/08/2017
% CALCULATE THE NUMERICAL MODEL FOR bigSCAle
% INPUT
% total_data: expression matrix
% edges: binning of expression values, internal variable
% ds: downsamples, can be 0 or an integer number. If ds=0 all cells are used to compute the model (small datasets of max 5k/10k cells).
% Otherwise ds is euqal to the number of cells used. 
% OUTPUT
% N: Raw enumerated counts for each combination
% N_pct: Enumerated counts transformed in cumulative distribution function

verbose=0;

% Performing down-sampling if ds>0
if (ds>0)
    num_samples=length(total_data(1,:))
    disp('Performing downsampling, as you requested');
    selected=randi(num_samples,ds,1);
    if issparse(total_data)
        total_data=full(total_data(:,selected));        
    else
        total_data=total_data(:,selected);
    end
   
end

% Remove genes with 1 or 0 umis/reads in the whole dataset. 
horrible =find(sum(total_data')<=1);

sprintf('I remove %g genes not expressed enough', length(horrible))
total_data(horrible,:)=[];


num_genes=length(total_data(:,1));
num_samples=length(total_data(1,:));

% normalize expression data for library size
somma_ex = sum(total_data);
total_data_lib_size=total_data./repmat(somma_ex,num_genes,1)* mean(somma_ex);

% old code, for debugging
%[~,~, central ]= population_calling( total_data_lib_size , 90);
%[ matrix ] = SC_log_transform( central, total_data_lib_size , 2 );
[ matrix ] = SC_log_transform( [], total_data_lib_size , 2 );
D = pdist(matrix','correlation');
Z = linkage(D,'ward');


% Adjusting max_group_size according to cell number
if (num_samples<1250) max_group_size=0.2; end
if (num_samples>=1250 && num_samples<=2000) max_group_size=0.17; end
if (num_samples>=2000 && num_samples<=3000) max_group_size=0.13; end  
if (num_samples>=3000 && num_samples<=5000) max_group_size=0.11; end
if (num_samples>=5000) max_group_size=0.9; end
    



disp('Calculating optimal cut of dendrogram for pre-clustering')
MAX_CUT_LEVEL=0.5;
cut_level=MAX_CUT_LEVEL; 
% Calculating optimal cutting depth for the pre-clustering
while 1
    while 1
        if (cut_level<0) break; end
        T = cluster(Z,'cutoff',cut_level*max(Z(:,3)),'criterion','distance');
        [counts]=hist(T,unique(T));
        if max(counts)<(max_group_size*num_samples) && length(unique(T))<(0.25*num_samples)
            break;
        end
        cut_level=cut_level-0.01;
    end
    if (cut_level>0) break;
    else
        max_group_size=max_group_size+0.01;
        cut_level=MAX_CUT_LEVEL;
    end
end

sprintf('Optimal cut found at %g %% with max_group_size %g %%',cut_level*100,max_group_size*100)
dendrogram(Z,Inf,'ColorThreshold', cut_level*max(Z(:,3)))
saveas(gcf,'./results/Pre_clustering.png')
close all

h = waitbar(0,'Calculating numerical model: ');

N_tot=zeros(length(edges),length(edges) );

% Enumerate the combinations
tot_computed_cells=0;
for k=1:max(T)
    waitbar(tot_computed_cells/num_samples,h,sprintf('Calculating numerical model: pool of %g cells',nnz(find(T==k))) );
    N = SC_1vs1_background( total_data(:,find(T==k)) , [] , edges );
    N_tot=N_tot+N;
    if verbose
        sprintf('Group %g, total of %g values',k,sum(sum(N_tot)))
    end
    tot_computed_cells=tot_computed_cells+nnz(find(T==k));
    waitbar(tot_computed_cells/num_samples);
end
close(h);
N=N_tot;
clear N_tot

% Calculating percentiles for the cumulative distribution function
for k=1:length(N)
    for j=1:length(N)
        N_pct(k,j) = sum(N(k,j:end))/sum(N(k,:));
    end
end
N_pct(isnan(N_pct))=1;

% old code for debugging
dereg_5=N_pct<0.05;
dereg_5=dereg_5+dereg_5';

dereg_1=N_pct<0.01;
dereg_1=dereg_1+dereg_1';


end

