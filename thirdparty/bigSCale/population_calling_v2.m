function [driving_genes genes_related score D]= population_calling_v2( total_data,option_coreg,min_driving,min_coreg)
% population_calling_v2 ************************************************************************
% GIOVANNI IACONO, CNAG, 16/08/2017
% Calculates overdispered genes which will be used to calculate cell-cell
% distances

% INPUT
% total_data: expression matrix
% bad_indexes: optional, force these genes to be excluded from
% overdispersed ones
% option_coreg: if 0, return overdispersed genes, if 1, selects only the
% overdispersed genes with at least min_coreg correlated genes.

% OUTPUT
% D_scores: Are the pvalues (uncorrected) in form of Z-scores 


score=[];
genes_related={};
num_genes=length(total_data(:,1))
num_samples=length(total_data(1,:))


% % Min Z-score to classify a gene as overdispersed
% min_driving=4 % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% %min_driving=1.645 % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% 
% 
% % Minimum number of co-regulated genes 
% % min_coreg=9 % linarsson
% min_coreg=5;



min_cells=max( [15 round(0.002*length(total_data(1,:)))] ) % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
skwed_cells=max([5 round(0.002*length(total_data(1,:))) ]) % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


% normalizing for library size
total_data=total_data ./ repmat(sum(total_data),length(total_data),1) * mean(sum(total_data));


% Removing too skewed genes
temp=sort(total_data,2);
A=log2(mean(temp(:,end-skwed_cells:end)'));
B=std(temp(:,end-skwed_cells:end)')./(mean(temp(:,end-skwed_cells:end)'));
[c ix]=sort(A);
blocco=1000;
conta=1;
for k=1:(length(A)-blocco)
    x_coord(conta)=mean(A(ix(k:k+blocco)));
	y_coord(conta)=mean(B(ix(k:k+blocco)));
	conta=conta+1;
end

okay=find(~isnan(x_coord) & ~isnan(y_coord));
f = fit(x_coord(okay)',y_coord(okay)','smoothingspline','SmoothingParam',0.25); %0.75

figure(1)
scatter(A,B,'.'); 
hold on
plot(f,x_coord',y_coord','k');
y = feval(f,A);
scatter(A,y,'.','r')
skewed_genes=find((B'./y)>4 & A'<prctile(A',99.5))
hold on
scatter(A(skewed_genes),B(skewed_genes),'r'); 
%scatter(A(6835),B(6835),'r','filled'); 
saveas(gcf,'./results/OD_genes_1_skewed.png')
clear A B f y x_* y_* conta
close all




okay=find( sum(total_data'>0)>min_cells );
sprintf('Using %g genes detected in at least >%g cells',length(okay),min_cells)

okay=setdiff(okay,skewed_genes);
sprintf('Further reducing to %g geni after discarding skewed genes', length(okay))

A=log2( mean(total_data(okay,:)') );
B=std(total_data(okay,:)')./mean(total_data(okay,:)');
  
% old code for debugging
% A=log2( mean(total_data') );
% B=std(total_data')./mean(total_data');


[c ix]=sort(A);
blocco=100;
conta=1;
for k=1:(length(A)-blocco)
    x_coord(conta)=mean(A(ix(k:k+blocco)));
	y_coord(conta)=mean(B(ix(k:k+blocco)));
	conta=conta+1;
end

f = fit(x_coord',y_coord','smoothingspline','SmoothingParam',0.75); %0.75

figure(2)
scatter(A,B,'.'); 
hold on
plot(f,x_coord',y_coord','k');
y = feval(f,A);
scatter(A,y,'.','r')
saveas(gcf,'./results/OD_genes_2_std.png')
close all;
% Part 2

clear x_coord y_coord
B_corr=B-y';
[c ix]=sort(A);
blocco=100;
conta=1;
for k=1:(length(A)-blocco)
    x_coord(conta)=mean(A(ix(k:k+blocco)));
	y_coord(conta)=trimm_std(B_corr(ix(k:k+blocco)));
	conta=conta+1;
end

f = fit(x_coord',y_coord','smoothingspline','SmoothingParam',0.75); %0.75
figure(3)
scatter(A,B_corr,'.'); 
hold on
plot(f,x_coord',y_coord','k');
y = feval(f,A);
scatter(A,y,'.','r')
saveas(gcf,'./results/OD_genes_3_mean.png')
close all;

score=zeros(num_genes,1); 
exprs=zeros(num_genes,1); 
score(okay)=B_corr./y';
exprs(okay)=A;
driving_genes=find(score>min_driving);

figure(4) % score e driving genes are indexed on the whole matrix
scatter(exprs,score,'.');
hold on
scatter(exprs(driving_genes),score(driving_genes),'.','r');
%scatter(exprs(686),score(686),'k','filled');


sprintf('Found %g overdispersed genes with Zscore>%.2f',length(driving_genes),min_driving)
% if any(bad_indexes)
%     sprintf('From which I remove %g bad indexes, remaining with %g driving genes',length(intersect(bad_indexes,driving_genes)), length(driving_genes)-length(intersect(bad_indexes,driving_genes)))
%     driving_genes=setdiff(driving_genes,bad_indexes);
% end


if option_coreg==0
    return;
end

alarm_big=0;
if num_samples>30000
    disp('Detecting big dataset, calculating correlation only for driving genes');
    okay=driving_genes;
    alarm_big=1;
end


matrix  = SC_log_transform( [], total_data(okay,:) , 2 );
disp('Calculating correlation ....');

D = squareform(pdist(matrix,'correlation'));
% Esimating correlation treshold using 10000 random correlations and
% selecting top 1% correlations
%p=randperm(numel(D));
%max_dist=prctile(nonzeros(D(p(1:10000))),0.025); %0.025
sprintf('Terminated calculating correlation')
%D = sqdist(matrix',matrix');

% calculating genes_related and tot_coreg
[sorted ix]=sort(D);

max_dist=0.9;
if isrow(okay)
    okay=okay';
end
conta=1;


while 1
    genes_related={};
    for k=1:length(sorted)
        genes_related{okay(k),1}=okay ( ix(find(sorted(:,k)>0 & sorted(:,k)<max_dist) , k) );
    end
    
    tot_coreg=cellfun(@length,genes_related);
    ok_coreg=find(tot_coreg>=min_coreg);
    sprintf('%g / %g driving genes present at leat %g genes co-regulated (max_dist=%.2f)',length(intersect(driving_genes,ok_coreg)),length(driving_genes),min_coreg,max_dist);
    
    % Taking all the genes related to the driving genes with tot_coreg>=min_coreg
    dummy = unique(cell2mat(genes_related(intersect(driving_genes,ok_coreg))));
    % To which I add  the same driving genes with tot_coreg>=min_coreg
    dummy = unique([dummy ; intersect(driving_genes,ok_coreg)]);
    sprintf('Found %g driving genes coreg with max_dist=%.2f',length(dummy),max_dist)
    max_dist=max_dist-conta/100
    if length(dummy)<=length(driving_genes)*1.2
        disp('Converged')
        break;
    end
end

% Restricting to the good driving genes
driving_genes=intersect(driving_genes,ok_coreg);

% Taking all the blocks of co-regulated genes associated with indexes of driving genes
dummy = unique(cell2mat(genes_related(driving_genes)));
if isrow(dummy)
    dummy=dummy';
end

% to which I add the same good driving genes
driving_genes= unique([dummy ; driving_genes]);
%driving_genes=setdiff(driving_genes,bad_indexes);
title(sprintf('Final number: %g driving genes',length(driving_genes)))
saveas(gcf,'./results/OD_genes_4_final.png')
close all;

sprintf('Final number: %g driving genes',length(driving_genes))

end


