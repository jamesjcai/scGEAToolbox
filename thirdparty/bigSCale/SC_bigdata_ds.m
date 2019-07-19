function [distances cells max_factor] = SC_bigdata_ds( total_data, N_pct , edges,somma_ex)
%
num_samples=length(total_data(:,1))

if num_samples<=100000
   SIZE=10000;
end
if num_samples>100000 &  num_samples<=300000
   SIZE=7500;
end
if num_samples>300000
   SIZE=5000;
end

max_factor=max(1,round(num_samples/20000));
max_factor=min(max_factor,8)

t=[];

num_genes=length(total_data(:,1))
num_samples=length(total_data(1,:))

% configuring dereg
soglia_dereg=0.01;

% symmetrizing 
[r c]=size(N_pct);
for k=1:r
    for h=1:c
        if k>h
            N_pct(k,h)=N_pct(h,k);
        end
    end
end

%  dereg
dereg=N_pct<soglia_dereg;


%  log_scores and infinites
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

% negative signes and zeroing diagonal
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
    
% cose richieste dal MEX log_scores e rimuovo infiniti
vector=([0:75000*10]); % puoi alzarlo se il codice si blocca
[~,ind_A]=histc(vector/10,edges);
ind_A=uint8(ind_A-1);
log_scores=abs(log_scores.*dereg); 
clear vector
%total_data_norm=total_data./repmat(somma_ex,num_genes,1);
    
conta=1;
link_wb=waitbar(0);


cells=zeros(num_samples,SIZE,'uint32'); %num_samples
distances=zeros(num_samples,SIZE,'single'); %num_samples


log_scores=single(log_scores);
somma_ex=single(somma_ex);

for h=1:num_samples
        
	waitbar(h / num_samples,link_wb,sprintf('%g / %g that is %.1f%%',h,num_samples,h / num_samples*100))
        
    to_test=setdiff( randi(num_samples,round(SIZE*1.5),1) , h ) ; % setdiff works also as unique
    to_test=to_test(randperm(length(to_test)));
    cells(h,:)=to_test(1:SIZE);
    
    total_data_norm=single(total_data(:,[h cells(h,:)]))./repmat(somma_ex([h cells(h,:)]),num_genes,1);
    
    distances(h,:)=SC_dist_MEX(total_data_norm,log_scores,ind_A,somma_ex,single(1)); %distances(h,:)=
end
    


end

