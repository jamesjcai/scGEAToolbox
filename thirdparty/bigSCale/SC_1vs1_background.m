function [ N N_pct ]=SC_1vs1_background( total_data, indexes , edges )
% SC_1vs1_background ************************************************************************
% GIOVANNI IACONO, CNAG, 16/08/2017
% ENUMERATES THE COMBINATIONS FOR bigSCale´s numerical model
% INPUT
% total_data: expression matrix
% edges: binning of expression values, internal variable
% indici: Default empty, use all genes. Otherwise it specifies the genes to be used.
% OUTPUT
% N: Raw enumerated counts for each combination
% N_pct: Enumerated counts transformed in cumulative distribution function

num_genes=length(total_data(:,1));
num_samples=length(total_data(1,:));
verbose=0;

if isempty(indexes)
    if verbose
        disp('Using all genes')
    end
    indexes=[1:num_genes];
end

% % normalize expression data for library size
somma_ex = sum(total_data);
total_data_norm=total_data./repmat(somma_ex,num_genes,1);


minimo_cellule=round(num_samples/40);

% Removing genes with low read counts ..
okay=find(sum(total_data'>0)>minimo_cellule);
okay=intersect(okay,indexes);
if verbose
sprintf('I reduce the analysis to %g genes expressed in at least %g cells',length(okay),minimo_cellule)
end

total_data_norm = total_data_norm(okay,:);

% Initiating variable N
N=zeros(length(edges));
if verbose
sprintf('Calculating %g couples',(num_samples*num_samples - num_samples)/2)
end

% Enumerating frequency of combinations by use of hist3
for i=1:num_samples
    for j=i:num_samples
        %idx=randi(num_samples,1,2);
        if ~(i==j)
            A= total_data_norm(:,i) * mean ( somma_ex( [i j] ) );
            B= total_data_norm(:,j) * mean ( somma_ex( [i j] ) );
            n=hist3([A B] ,'Edges' ,  {edges edges});
            N=N+n;
        end
    end
end
N=N+N';


% Calculating percentiles
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














% old data for debug

% for k=1:length(N)
%     dummy=min(find(N_pct(k,:)<pval));
%     if any(dummy)
%         position(k)=min(find(N_pct(k,:)<pval));
%     else
%         position(k)=NaN;
%     end
% end



% for j=1:1
%     j
%     idx=randi(1847,1,2);
% 
%     A=log2( total_data_norm(:,idx(1)) * mean ( somma_ex( [idx(1) idx(2)] ) ) );
%     B=log2( total_data_norm(:,idx(2)) * mean ( somma_ex( [idx(1) idx(2)] ) ) );
%     okay=~isinf(A) & ~isinf(B);
%     A=A(okay);
%     B=B(okay);
%     [x_rotated y_rotated] = rotate(A',B');
%     [x_rotated ix]=sort(x_rotated);
%     y_rotated = y_rotated(ix);
%     clear ix
%     block=50;
%     clear fitting pos
%     for k=1:length(x_rotated)-block
%         fitting(k)=mean(abs(y_rotated(k:k+block) ) );
%         pos(k)=mean( x_rotated(k:k+block) );
%     end
%     plot(pos,fitting,'r');
%     hold on
%     %scatter(x_rotated,abs(y_rotated))
% end



% % figure(1)
% % hold on
% 
% %edges=min(X(~isinf(X))): (max(X(~isinf(X)))-min(X(~isinf(X))))/10 :max(X(~isinf(X)));
% edges=[-3: 13/100 :10];
% 
% conta=1;
% 
% for j=1:1000
%     j
%     idx=randi(1847,1,2);
% 
%     A=total_data_norm(:,idx(1));
%     B=total_data_norm(:,idx(2));
% 
%     X = log2(  mean([A B]') * mean ( somma_ex( [idx(1) idx(2)] ) ) );
%     Y =  abs((A-B)')./mean([A B]');
%     %Y =  std([A B]')./mean([A B]');
%     
%     %okay=find(~isinf(X) & ~isnan(X));
%     %tot_X(conta:conta+length(okay)-1)=X(okay);
%     %tot_Y(conta:conta+length(okay)-1)=Y(okay);
%     %conta=conta+length(okay);
%     for k=1:length(edges)-1
%         fitting(j,k)=nanmean(Y ( X>=edges(k) & X<edges(k+1) ));
%     end
%     
% end


    
% end

