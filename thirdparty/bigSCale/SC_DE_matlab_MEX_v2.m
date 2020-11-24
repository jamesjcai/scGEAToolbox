function [ D_scores F_changes D_scores_WC D_scores_fusion]=SC_DE_matlab_MEX_v2( total_data, N_pct , edges,indexes, name_file_xls,kg_or_ens,samples_names,speed_preset)
% SC_DE_matlab_MEX_v2 ************************************************************************
% GIOVANNI IACONO, CNAG, 16/08/2017
% DIFFERENTIAL EXPRESSION WITH bigSCale
% Compares two groups of cells to find differentially expressed genes
% INPUT
% total_data: expression matrix
% N_pct: numerical model variable
% edges: binning of expression values
% indices: indexes of the two groups to compare
% nome_file_xls: optional, to write output to excel file
% kg_or_ens: optional, to write output to excel file
% samples_names: optional, to write output to excel file
% OUTPUT
% D_scores: Are the pvalues (uncorrected) in form of Z-scores 

% MARIONI_norm=0;
MARIONI_norm=0;

speed_preset_found=0;

if isequal(speed_preset,'sfast')
    max_size=1500;
    speed_preset_found=1;
end
if isequal(speed_preset,'fast')
    max_size=1500;
    speed_preset_found=1;
end
if isequal(speed_preset,'normal')
    max_size=3000;
    speed_preset_found=1;
end
if isequal(speed_preset,'slow')
    max_size=5000;
    speed_preset_found=1;
end

if speed_preset_found==0
   error('Wrong speed preset'); 
end

% Preliminary cleaning for very large size comparisons *******************
if length(indexes{1})>max_size
    temp=indexes{1};
    temp=temp(randperm(numel(temp)));
	indexes{1}=temp(1:max_size);
end
if length(indexes{2})>max_size
    temp=indexes{2};
    temp=temp(randperm(numel(temp)));
	indexes{2}=temp(1:max_size);
end
% ************************************************************************


num_genes=length(total_data(:,1))
num_samples=length(total_data(1,:))


total_data_safe=total_data;

% remove the genes with zero reads
geni_okay = find(sum(total_data(:,[indexes{1} ; indexes{2}])')>0 );
total_data=total_data(geni_okay,:);

% 1) MAKE N_pct symmetric
[r c]=size(N_pct);
for k=1:r
    for h=1:c
        if k>h
            N_pct(k,h)=N_pct(h,k);
        end
    end
end

% 3) Calculate log_scores
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

% 4) Put minus signs and set diagonal to zero
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
    
sum_ex=sum(total_data);
if MARIONI_norm==0
    total_data_norm=total_data./repmat(sum_ex,length(total_data(:,1)),1);
else
    total_data_norm = SC_normalizza_scran( total_data );
    total_data_norm=total_data_norm./mean(sum(total_data_norm)); % Otherwise crushes C, total_data_norm here is supposed to be a fraction !
end





vector=([0:1000000*10]); % Increase this number if there are errors when runing the code
[~,ind_A]=histc(vector/10,edges);
ind_A=uint8(ind_A-1);
%tic; 
dim=mean(cellfun(@length,indexes))

if isequal(speed_preset,'sfast')
     max_rep =  1
end

if isequal(speed_preset,'fast')
    if dim < 200
     max_rep =  5 
    end
    if dim >=200 & dim < 500
     max_rep =  3  
    end
    if dim >= 500 & dim<1000
     max_rep =  2 
    end
    if dim>=1000
     max_rep =  1
    end
end

if isequal(speed_preset,'normal')
    if dim < 500
     max_rep =  5 
    end
    if dim >=500 & dim < 1000
     max_rep =  3  
    end
    if dim >= 1000 & dim<2000
     max_rep =  2 
    end
    if dim>=3000 
     max_rep =  2
    end
end

if isequal(speed_preset,'slow')
    if dim < 500
     max_rep =  5 
    end
    if dim >=500 & dim < 1000
     max_rep =  5  
    end
    if dim >= 1000 & dim<2000
     max_rep =  4 
    end
    if dim>=2000 & dim<4000
     max_rep =  2
    end
    if dim>=4000
     max_rep =  1
    end
end



F_changes=log2(mean(total_data_norm(:,indexes{1})')./mean(total_data_norm(:,indexes{2})'));
dummy=zeros(num_genes,1);
dummy(geni_okay)=F_changes;
F_changes=dummy;
clear dummy

disp('Calculating scores');
for repetitions=1:max_rep
    % calculating reshuffled DEs
	repetitions
	idx=union(indexes{1},indexes{2});
	idx=idx(randperm(numel(idx)));
	fake_A=idx(1:length(indexes{1}));
	fake_B=idx(length(indexes{1})+1:end);  
	[D_scores_fake(:,repetitions) comp_num_fake(:,repetitions)]=SC_DE_MEX_double(total_data_norm(:,fake_A),log_scores,ind_A,sum_ex(fake_A),total_data_norm(:,fake_B),sum_ex(fake_B));  
end

D_scores_fake=reshape(D_scores_fake,numel(D_scores_fake),1);
comp_num_fake=reshape(comp_num_fake,numel(comp_num_fake),1);

% SC_DE_MEX_double() Mex C Matlab CODE ***********************
% Here the C code embebbed into Matlab starts to perform the
% computationally intense part
[D_scores comp_num]=SC_DE_MEX_double(total_data_norm(:,indexes{1}),log_scores,ind_A,sum_ex(indexes{1}),total_data_norm(:,indexes{2}),sum_ex(indexes{2})); 
%t=toc

dummy=zeros(num_genes,1);
dummy(geni_okay)=D_scores;
D_scores=dummy;
clear dummy
dummy=zeros(num_genes,1);
dummy(geni_okay)=comp_num;
comp_num=dummy;
clear dummy


% Calculate x_coor e y_coord
A=log2(comp_num_fake+1);
B=D_scores_fake;

[c ix]=sort(A);
%scatter(A,B,'b','filled','MarkerFaceAlpha',0.1)
blocco=100;
conta=1;
for k=1:round(length(c)/200):length(c)-blocco
	x_coord(conta)=mean(c(k:k+blocco));
	y_coord(conta)=std(B(ix(k:k+blocco)));
    conta=conta+1;
end
[ix]=find(~isnan(y_coord));

% fitting eval
f = fit(x_coord(ix)',y_coord(ix)','smoothingspline','SmoothingParam',0.97);

A=log2(comp_num+1);

yy = feval(f,A);

% DETERMINE tresholds
tresh1=log2((length(indexes{1})*length(indexes{2}))/6+1) % 6=16.7%
min_cells=15;
tresh2=max([  log2( (min_cells*length(indexes{1})) +1 )   log2( (min_cells*length(indexes{2})) +1 ) ] )
tresh=max(tresh1,tresh2)
    
% CORRECT D_scores
[~,pos]=min(abs(A-tresh));
yy(A<=tresh)=mean( yy(pos) );

% old code, for debugging
%hold on;
%plot(x_coord,y_coord,'k');
%scatter(A,yy,'r','.');
%figure(2)
%scatter(A,D_scores);

D_scores=D_scores./yy;


% RUN Wilcoxon


D_scores_WC=zeros(length(geni_okay),1);
for k=1:length(geni_okay)  
    D_scores_WC(k)=ranksum(total_data_norm(k,indexes{1}), total_data_norm(k,indexes{2}));
end

D_scores_WC=D_scores_WC/2;
D_scores_WC(isnan(D_scores_WC))=1;
D_scores_WC=-pval_to_Zscore(abs(D_scores_WC),1);

dummy=zeros(num_genes,1);
dummy(geni_okay)=D_scores_WC;
D_scores_WC=dummy;
D_scores_WC=D_scores_WC.*sign(D_scores);
clear dummy


factor=max(D_scores_WC)/max(D_scores);
sprintf('Normalizing WC over BS Zscores: factor %.1f (max(BS)=%.1f, max(WC)=%.1f)',factor,max(D_scores),max(D_scores_WC))
D_scores_WC=D_scores_WC/factor;

sprintf('after  max(WC)=%.1f',max(D_scores_WC))

D_scores_fusion=sqrt(D_scores.^2+D_scores_WC.^2).*sign(D_scores);

D_scores=D_scores_fusion;


% old code, for debugging
%figure(3)
%scatter(A,D_scores);
%hold on
%catter(A(2691),D_scores(2691),'y','filled');
%sprintf('Ci ho messo %g minuti (ovvero %g secondi)',round(t/60),round(t) )






if  length(name_file_xls)
    
    total_data=total_data_safe;
    total_data_norm=total_data./repmat(sum_ex,num_genes,1)*mean(sum_ex);
    f1=sum(total_data_norm(:,indexes{1})'>0)./sum(total_data_norm(:,indexes{1})'>=0);
    f2=sum(total_data_norm(:,indexes{2})'>0)./sum(total_data_norm(:,indexes{2})'>=0);
    
    e1=zeros(size(f1));
    e2=zeros(size(f1));
    
    disp('Preparing data for exporting ....')
    for k=1:num_genes
        e1(k)=mean(nonzeros(total_data_norm(k,indexes{1})));
        e2(k)=mean(nonzeros(total_data_norm(k,indexes{2})));
        
        if isnan(e1(k))
            e1(k)=0;
        end
         if isnan(e2(k))
            e2(k)=0;
         end
        
        if f1(k)>f2(k)
            e2(k) = (f2(k)*e2(k)) / (f1(k)) ;
        else
            e1(k) = (f1(k)*e1(k)) / (f2(k)) ;
        end
    end
    
    
    disp('Exporting ...... ')
    header={ 'GENE ID'	'GENE NAME'	'GENE TYPE' 'DESCRIPTION' 'D_score'	sprintf('Freq.%s',samples_names{1})	 sprintf('Exprs.%s',samples_names{1})	sprintf('Freq.%s',samples_names{2})	 sprintf('Exprs.%s',samples_names{2})	sprintf('Log2(%s / %s)',samples_names{2},samples_names{1}) }
    [ ~, ix]=sort(D_scores);
    to_write=[ kg_or_ens(:,1:end) mio_mat2cell([ D_scores f1' e1' f2' e2' log2(e2'./e1')])];    
    delete(sprintf('./../data/SC_output_compare_gene_expr/%s.xlsx',name_file_xls))
    copyfile('./../data/SC_output_compare_gene_expr/TEMPLATE.xlsx',sprintf('./../data/SC_output_compare_gene_expr/%s.xlsx',name_file_xls))
    xlswrite(sprintf('./../data/SC_output_compare_gene_expr/%s.xlsx',name_file_xls),header,1,'A1');
    xlswrite(sprintf('./../data/SC_output_compare_gene_expr/%s.xlsx',name_file_xls),to_write(ix,:),1,'A2');    
    
end


end









