function [ D_ls ]=SC_1vs1_MEX( total_data, N_pct , edges,sum_ex , remove_genes)
% SC_1vs1_MEX ************************************************************************
% GIOVANNI IACONO, CNAG, 16/08/2017
% DISTANCE MATRIX
% Calculates the distance matrix
% INPUT
% total_data: expression matrix
% N_pct: numerical model variable
% edges: binning of expression values
% indices: indexes of the two groups to compare
% sum_ex: library size
% remove genes: genes of the confunding signatures which effect must be remved or reduced. These genes are obtained via the fuction SC_bool_v2. 
% OUTPUT
% D_ls: Distance matrix


num_genes=length(total_data(:,1))
num_samples=length(total_data(1,:))

 if issparse(total_data)
    total_data=full(total_data);
 end

% Configure deregulated treshold 
soglia_dereg=0.01;

% 1) Making N_pct symmetric
[r c]=size(N_pct);
for k=1:r
    for h=1:c
        if k>h
            N_pct(k,h)=N_pct(h,k);
        end
    end
end

% 2) Finding deregulated
dereg=N_pct<soglia_dereg;

% 3) Calculating log_scores and removing infinite values
log_scores=-log10(N_pct);
[r c]=size(log_scores);
for k=1:r
    massimo=max(log_scores(k,~isinf(log_scores(k,:))));
    for h=1:c
        if isinf(log_scores(k,h))
            log_scores(k,h)=massimo;
        end
    end
end

% 4) Putting also negative signs and zeroing the diagonal 
    for k=1:r
        for h=1:c
               if k>h 
                   log_scores(k,h)=-log_scores(k,h);
               end
               if k==h 
                   log_scores(k,h)=0;               
               end
        end
    end



total_data_norm=total_data./repmat(sum_ex,num_genes,1);
    
vector=([0:1000000*10]); % increase if code blocks
[~,ind_A]=histc(vector/10,edges);
ind_A=uint8(ind_A-1);
log_scores=abs(log_scores.*dereg); 

if any(remove_genes)
    disp('Detect that you want to remove genes')
    disp('NOTE THAT GENE IDX MUST BE MAPPED OVER THE OVERDISPERSED, FIX THIS!')
    pause(5);
    vector_weights=ones(num_genes,1)';
    vector_weights(remove_genes)=0;
    tic; D_ls=SC_dist_MEX_double_rv(total_data_norm,log_scores,ind_A,sum_ex,0,vector_weights); t=toc
else
    tic; D_ls=SC_dist_MEX_double(total_data_norm,log_scores,ind_A,sum_ex,0); t=toc
end

if iscolumn(D_ls)
    D_ls=D_ls';
end

sprintf('Time elapsed %g minutes (that is %g seconds)',round(t/60),round(t) )

end









