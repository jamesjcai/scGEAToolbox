
clear;  clc;
tic

data_dir='data/'%'D:/AanchalMongia_phdClg/Phd/DATA_SINGLE_CELL/'
addpath(genpath('Dependencies')); 
addpath(genpath('functions')); 

dataname='Preimplantation'
  
 %% Data read
if(strcmp(dataname, 'Usoskin') | strcmp(dataname,'Zeisel') )
    data=csvread([data_dir 'raw_data/' dataname '_raw_data.csv'],1,1)' ;
else
   sample_dir = [data_dir 'raw_data/' dataname '_dataset/hg19'];
   [data, gene_names, gene_ids, cells] = read_raw_10x( strcmp(dataname,'Preimplantation') ,sample_dir);% 0 for jurkat,1 for zygote
end

%% CALL MCIMPUTE
[data_recovered,~,data_recovered_raw, gene_names]=call_mcImpute(data,'dataname',dataname);%,'gene_names',gene_names,'gene_ids', gene_ids,'pro_dir',pro_dir);

%% Save results
mkdir(['RecoveredMatrices/' dataname '_imputed']); 
save(['RecoveredMatrices/' dataname '_imputed/rec.mat'],'data_recovered') 
time_taken=toc 



