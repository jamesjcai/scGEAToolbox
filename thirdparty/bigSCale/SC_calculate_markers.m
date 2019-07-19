function [markers_hard markers_soft scores_hard scores_soft I_scores FC_hard FC_soft]= SC_calculate_markers( total_data, indices, N_pct,treshold,edges,speed_preset)
% SC_calculate_markers ************************************************************************
% GIOVANNI IACONO, CNAG, 16/08/2017
% calculates scores for markers, internal function

% INPUT
% total_data: expression matrix, not normalized
% indices: indices of the clusters
% N_pct: Model variable
% treshold: internal variable
% edges: internal variable

% OUTPUT
% I_scores: internal variable needed by SC_exprs_boolean to write to disk
% hierarchical markers

sum_ex=sum(total_data);

sprintf('I use speed_preset = %s',speed_preset)


tic;
h=waitbar(0,'Calculating markers: 0.0%');
conta=0;
to_calculate=length(indices)*(length(indices)-1)/2

if ~iscell(N_pct)

    for j=1:length(indices)-1 % run trough all the pairwise comparision between clusters
        for k=j+1:length(indices)
            tic;
            sprintf('j=%g, k=%g, %g vs %g cells',j,k,numel(indices{j}),numel(indices{k}))
            % For each comparison perform DE and save the data in I_scores
            I_scores{j,k}=SC_DE_matlab_MEX_v2( total_data, N_pct , edges, indices([j k]),[],[],[],speed_preset ); % SC_1vs1( total_data,indici([j k]), N_pct, somma_ex,[],[],[]);
            I_scores{k,j}=-I_scores{j,k};

            F_change{j,k}=SC_DE_fc(total_data,indices([j k]));
            F_change{k,j}=-F_change{j,k};
            
            conta=conta+1;
            waitbar(conta/to_calculate,h,sprintf('Calculating markers ... %.1f %%',conta/to_calculate*100));
            
            t=toc
        end
    end

close(h);
 
else
    % Calculate fold_changes too... 
    I_scores=N_pct;
    for j=1:length(indices)-1
        for k=j+1:length(indices)  
            sprintf('j=%g, k=%g',j,k)
            F_change{j,k}=SC_DE_fc(total_data,indices([j k]));
            F_change{k,j}=-F_change{j,k};   
        end
    end
end

% cycle all the clusters of cell
for j=1:length(indices)
    
    % calculate block_values
    block_values=[];
    for k=1:length(indices)
        if ~(k==j)    
            block_values=[block_values I_scores{j,k} ];
        end
    end
    block_values=-block_values;
    
   % block_values il block_fc
    block_fc=[];
    for k=1:length(indices)
        if ~(k==j)    
            block_fc=[block_fc F_change{j,k} ];
        end
    end
    block_fc=-block_fc;
    
    % calculate lv.1 markers (markers hard)

    ix_good = find(sum(block_values'>treshold)==(length(indices)-1));
    score=min(block_values(ix_good,:)');%median(blocco_valori(ix_good,:)');
    [c ix]=sort(score,'descend');
    markers_hard{j}= ix_good(ix); 
    scores_hard{j}= score(ix)';
    
    score=min(block_fc(ix_good,:)');
    FC_hard{j}= score(ix)';

    if length(indices)>3
    % calculate lv.2 markers (markers soft)
    [blocco_valori_sorted ix_bv] = sort(block_values');
    dummy_scores_soft=min(blocco_valori_sorted(2:end,:))';%median(blocco_valori_sorted(2:end,:))';
    ix_good = find(sum(block_values'>treshold)>=(length(indices)-2));
    score=dummy_scores_soft(ix_good);
    [c ix]=sort(score,'descend');
    markers_soft{j}= ix_good(ix); 
    scores_soft{j}= score(ix)';
    
    ix_bv=ix_bv';
    for k=1:length(ix_bv(:,1))
        block_fc(k,:)=block_fc(k,ix_bv(k,:));
    end
    score=min(block_fc(ix_good,2:end)');
    FC_soft{j}= score(ix)';
    else
        markers_soft=[];
        scores_soft=[];
        FC_soft=[];
    end
    
    
end

for k=1:length(markers_soft)
    for j=1:length(markers_soft)
        common{j,k}=intersect(markers_soft{j},markers_soft{k});
    end
end


total_time=toc;
sprintf('Calculate markers: Elapsed %g minutes ,%g seconds',round(total_time/60),round(total_time) )

end

