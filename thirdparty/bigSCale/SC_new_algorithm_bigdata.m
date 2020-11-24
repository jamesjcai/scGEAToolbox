function [N N_pct difference] = SC_new_algorithm_bigdata( total_data , edges)
% SC_new_algorithm_bigdata ************************************************************************
% GIOVANNI IACONO, CNAG, 16/08/2017
% CALCULATE THE NUMERICAL MODEL FOR bigSCAle for very big datasets.
% It iteratevely lunches SC_new_algorithm with random downsamples cells until convergence.
% INPUT
% total_data: expression matrix
% edges: binning of expression values, internal variable
% OUTPUT
% N: Raw enumerated counts for each combination
% N_pct: Enumerated counts transformed in cumulative distribution function
% difference: For debugging

N_pct=[];
conta=1;
margine=20;
stop=95; % when stop% of values has en error below margin% stop.
max_cicli=8;

while 1
    [N(:,:,conta) N_pct(:,:,conta) ] = SC_new_algorithm( total_data , edges, 5000 );

    if conta>1
        
        [ r c z]=size(N_pct);
        
        for i=1:r
            for j=1:c
                mean_value(i,j)=mean(N_pct(i,j,:));
                s_dev(i,j)=std(N_pct(i,j,:));
            end
        end
        
        sem=s_dev/sqrt(conta);
        difference=sem./mean_value*100;
        
        sprintf('Cycle %g, found %g/%g (%g %%) values within the limit of %g %%',conta,nnz(difference<margine),numel(difference),nnz(difference<margine)/numel(difference)*100,margine)
        pause(5);
        
    end

    

    if (conta>1 & nnz(difference<margine)/numel(difference)*100 > stop) | conta>=max_cicli
        break;    
    end
    
    conta=conta+1
    
end

disp('Calculating final values');
N=sum(N,3);
clear N_pct
% ricavo i percentili
for k=1:length(N)
    for j=1:length(N)
        N_pct(k,j) = sum(N(k,j:end))/sum(N(k,:));
    end
end
N_pct(isnan(N_pct))=1;

end
