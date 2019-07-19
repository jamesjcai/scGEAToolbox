function [all_markers matrix_raw_output T6_lists T12_lists T24_lists]=SC_bool_v2( I_scores , kg_or_ens, treshold,option)
% SC_bool_v2 ************************************************************************
% GIOVANNI IACONO, CNAG, 16/08/2017
% writes to disk an excel file with the hierarchal markers
% INPUT
% I_scores: output of SC_calcola_markers
% kg_or_ens: annotatation file with gene names
% SOGLIA,opzione: variable for SC_exprs_boolean

% OUTPUT
% all_markers: list of all significant markers
% matrix_raw_output: matrix containing the Zscores for all genes and all
% clusters. A high Zscore indicated that the given gene is likley to be a
% marker of the given cluster. Zscores can be changed into pvalues with the
% function Zscore_to_pval. The matrix can then be exported with dlmwrite.
% T6,T12,T24 lists: moodules of coexpressed genes (6 module, 12 modules of
% 24 modules). TO be used for batch effect removal.

% Maxium hierarchical level to look for. Default should be total clusters - 1
max_depth=length(I_scores)-1;
%max_depth=7;

% delete('./../data/SC_output_matrix_bool/matrix_bool.xlsx');
% copyfile('./../data/SC_output_matrix_bool/BOOL_model.xlsx','./../data/SC_output_matrix_bool/matrix_bool.xlsx')

delete('./results/matrix_bool.xlsx');
copyfile('BOOL_model.xlsx','./results/matrix_bool.xlsx')

used=[];
% at every cycle calls SC_exprs_boolean to calculate the markers of level "depth"

%max_depth=5;

for depth=1:max_depth
     depth
%     if depth<=1 
%         treshold_use=4;
%     else
%         treshold_use=treshold;
%     end
    [matrix_raw_output(:,:,depth) , ~, lists] = SC_exprs_boolean( I_scores , kg_or_ens, treshold,option,depth);
    used=[used ; cell2mat(lists')];
end
matrix_raw_output=mean(matrix_raw_output,3);
all_markers=unique(used);



% CALCULATE COEXPRESSED SIGNATURES USING I_SCORES
% too slow when lots of genes....
input=I_scores;
for i=1:length(input)
for j=1:length(input)
    if j>=i
        input{i,j}=[];
    end
end
end
input=cell2mat(reshape(input,1,numel(input)));

okay=find(sum(abs(input)>treshold,2));
input=input(okay,:);
disp('Calculating distance');
D = pdist(input,'cosine');


% %%CALCULATE COEXPRESSED SIGNATURES USING MATRIX_RAW_OUTPUT
% faster than the other ...
% okay=find(sum(matrix_raw_output,2));
% matrix_raw_output=matrix_raw_output(okay,:);
% disp('Calculating distance')
% D = pdist(matrix_raw_output,'cosine');



disp('Calculating dendrogram')
Z = linkage(D,'ward');

disp('Calculating lists')
T6 = cluster(Z,'MaxClust',6); 
T12 = cluster(Z,'MaxClust',12); 
T24 = cluster(Z,'MaxClust',24);

for k=1:max(T6)
T6_lists{k}=okay(find(T6==k));
end

for k=1:max(T12)
T12_lists{k}=okay(find(T12==k));
end

for k=1:max(T24)
T24_lists{k}=okay(find(T24==k));
end



disp('WRITING OUTPUT SIGNATURES')
delete('./results/output_signatures.xlsx');
copyfile('template_output_signatures.xlsx','./results/output_signatures.xlsx');
for k=1:length(T6_lists)
    name=xlsColNum2Str(k);
    xlswrite('./results/output_signatures.xlsx',sort(kg_or_ens(T6_lists{k})),'T6 lists',sprintf('%s2',name{1})); 
end
for k=1:length(T12_lists)
    name=xlsColNum2Str(k);
    xlswrite('./results/output_signatures.xlsx',sort(kg_or_ens(T12_lists{k})),'T12 lists',sprintf('%s2',name{1})); 
end
for k=1:length(T24_lists)
    name=xlsColNum2Str(k);
    xlswrite('./results/output_signatures.xlsx',sort(kg_or_ens(T24_lists{k})),'T24 lists',sprintf('%s2',name{1})); 
end




% % ALTERNATIVE OUTPUT
% 
% % average
% %output=Zscore_to_pval(mean(matrix_raw_output,3));
% 
% % specific level
% output=Zscore_to_pval(matrix_raw_output(:,:,5));
% 
% 
% 
% [r c]=size(kg_or_ens);
% if c==1
%     header={'GENE'};
% else
%     header={'GENE_ID' 'GENE NAME' 'DESCRIPTION'};
% end
% for j=1:length(I_scores)
%     header{1,end+1}=sprintf('C_%g',j);
% end
% xlswrite('./../data/SC_output_matrix_bool/matrix_bool.xlsx',header,'summary','A1');
% xlswrite('./../data/SC_output_matrix_bool/matrix_bool.xlsx',[kg_or_ens  mio_mat2cell(output)],'summary','A2');

end


