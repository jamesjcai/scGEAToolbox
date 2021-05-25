function [res_pos,res_neg,res_descr,p_gene] = MrGSEAPreranked(data_ranks,genelist,GeneSet,GeneSetName,file_name,opts)
% Computing Gene Set Enrichment Analysis (GSEA) with various ranking metrics.
% Major adjustments:
% - 12 gene ranking metrics implemented (plus its ranks and/or absolute values)
% - possibility to use external ranking metric method
% - possibility to use external Gene Set Database
% - phenotype permutation available for all ranking methods (even external)
% - increased efficiency by parallel computing
%
% Input:
% data - gene expression matrix [mxn]
% group - sample labels [1xn]
% data_entrez - gene entrez IDs [mx1]
% file_name - save file name 
% opts - structure of parameters
% Output:
% res_pos - table with positive GSEA results
% res_neg - table with negative GSEA results
% res_desc - description of tables with results
% p_gene - genes p-values from permutation test
%
% Author:
% Michal.Marczyk@polsl.pl

%check inputs
if nargin < 6
    opts = default_GSEA_opts();
end

% if nargin < 3
%     error('Too few input arguments.')
% end
% 
% if length(unique(data_entrez)) ~= length(data_entrez)
%     error('Non-unique gene identifiers found.')
% end

%calculate ranks
% [data_ranks,opts.tail] = calc_ranks(data,group,opts.rank_type,opts.abs,opts.tied_rank);

opts.tail='both';

%sort data by ranks
[data_ranks,ind] = sort(data_ranks,opts.sort_type);
genelist = genelist(ind);

% data = data(ind,:);
n = length(data_ranks);

% %load GS database
% if opts.show; disp(['Load ' opts.GS_name ' database.']); end
% tmp = load(opts.GS_name);
% GS = tmp.GS;

GS.entrez_nb=length(GeneSetName);
%filter database by removing genes in GS that are not present in data
ind_d = cell(1,length(GeneSetName));  %indices of all GS genes in dataset
RI = ind_d;
for a=1:length(GeneSetName)
    [common,ind_d{a}] = intersect_fast(genelist,string(GeneSet{a}));  %find genes that are in given GS
    %if length(common) ~= GS.entrez_nb(a)
        GS.entrez{a} = common;
        GS.entrez_nb(a) = length(common);    %find genes with ranks within GS
        GS.ID{a}=a;
        GS.descr{a}=GeneSetName{a};
    %end
    RI{a} = data_ranks(ind_d{a});
end
    
%filter database by min. number of genes in GS
if opts.GS_filt(1) > 0 || opts.GS_filt(2) < Inf
    del = GS.entrez_nb < opts.GS_filt(1) | GS.entrez_nb > opts.GS_filt(2);
    GS.ID(del) = [];
    GS.descr(del) = [];
    GS.nb = length(GS.ID);
    GS.entrez(del) = [];
    GS.entrez_nb(del) = [];
    ind_d(del) = [];
    RI(del) = [];
    if opts.show; disp([num2str(sum(del)) ' GS filtered, ' num2str(GS.nb) ' GS left in database.']); end
else
    if opts.show; disp([num2str(GS.nb) ' GS in database found.']); end
end

%pre-calculate ranks for permutation test
if opts.show
    disp('Calculate ranking metrics for all permutations.')
end
% ranks_perm = pre_calc_ranks(data_ranks,data,group,opts);

m=length(data_ranks);
ranks_perm = zeros(m,opts.perm_nb);
% rank_type = opts.rank_type;
% if_abs = opts.abs;
% if_trank = opts.tied_rank;
%if strcmp(opts.perm_type,'entrez')
    for a=1:opts.perm_nb
        ind = RandPermFast(m);
        ranks_perm(:,a) = data_ranks(ind);
    end

switch opts.tail
    case 'both'
        p_gene = 1-(sum(repmat(abs(data_ranks),1,opts.perm_nb) > ranks_perm,2))/opts.perm_nb;
    case 'left'
        p_gene = 1-(sum(repmat(data_ranks,1,opts.perm_nb) < ranks_perm,2))/opts.perm_nb;
    case 'right'
        p_gene = 1-sum(repmat(data_ranks,1,opts.perm_nb) > ranks_perm,2)/opts.perm_nb;
end
p_gene(ind) = max(p_gene,1/opts.perm_nb);

% separate calculations for each GeneSet
out = zeros(GS.nb,9); NES_perm = zeros(GS.nb,opts.perm_nb);
LEF = cell(GS.nb,1);
if opts.show
    [~,~] = mkdir('GSEA_plots');
%     fprintf('Progress:\n');
%     fprintf(['\n' repmat('.',1,GS.nb) '\n\n']);
end
parfor a=1:GS.nb 
    opts_tmp = opts;
    out_tmp = zeros(1,9);
    GS_tmp = GS;
        
    %sort ranks in GS to find LEF
    RI_tmp = sort(RI{a},opts_tmp.sort_type); 
    
    %Enrichment Score calculation
    Phit = zeros(n,1); 
    Pmiss = ones(n,1) * 1/(n-GS_tmp.entrez_nb(a));
    Phit(ind_d{a}) = (abs(RI{a}).^opts_tmp.p)/sum(abs(RI{a}));
    Pmiss(ind_d{a}) = 0;
    ES = cumsum(Phit - Pmiss);    
    [~,ind] = max(abs(ES));     %find maximum ES 
    ESmax = ES(ind);
    
% -- MOVE down --  
%     if opts.show
%         plot_classic(ES,ESmax,data_ranks,ind_d{a},GS.ID{a},Phit,Pmiss);
%     end
    
    % Find p-value by permutation test
    ES_perm = ES_perm_test(ranks_perm,genelist,GS_tmp.entrez{a},opts_tmp.perm_nb,opts_tmp.sort_type);
        
    out_tmp(1) = a;        % GS name  
    out_tmp(2) = GS_tmp.entrez_nb(a);     % GS Size  
    out_tmp(3) = ESmax;           % Enrichment Score 
    
    if ESmax>0    %positive ES
        [LEF_tmp,ind_LEF] = intersect_fast(GS_tmp.entrez{a},genelist(1:ind));   % genes from LEF 
        LEF_tmp(:,2) = RI_tmp(ind_LEF);      % ES from genes of LEF   
        ind = ES_perm > 0;
        tmp_NES = zeros(1,opts.perm_nb);
        tmp_NES(ind) = ES_perm(ind)/mean(ES_perm(ind));
        NES_perm(a,:) = tmp_NES;
        if ~isempty(ES_perm(ind))
            out_tmp(4) = max(1,sum(ES_perm >= ESmax))/length(ES_perm(ind));  % p-value
            out_tmp(5) = ESmax/mean(ES_perm(ind)); % Normalized Enrichment Score 
        else
            out_tmp(4) = NaN;  % p-value
            out_tmp(5) = NaN; % Normalized Enrichment Score 
        end
        
    else    %negative ES
        [LEF_tmp,ind_LEF] = intersect_fast(GS_tmp.entrez{a},genelist(ind:end));   % genes from LEF 
        LEF_tmp(:,2) = RI_tmp(ind_LEF);      % ES from genes of LEF   
        ind = ES_perm < 0;
        tmp_NES = zeros(1,opts.perm_nb);
        tmp_NES(ind) = ES_perm(ind)/abs(mean(ES_perm(ind)));
        NES_perm(a,:) = tmp_NES;
        if ~isempty(ES_perm(ind))
            out_tmp(4) = max(1,sum(ES_perm <= ESmax))/length(ES_perm(ind));  % p-value
            out_tmp(5) = ESmax/abs(mean(ES_perm(ind))); % Normalized Enrichment Score 
        else
            out_tmp(4) = NaN;  % p-value
            out_tmp(5) = NaN; % Normalized Enrichment Score 
        end
    end
    out_tmp(9) = size(LEF_tmp,1);          % observed counts in LEF
    
    out(a,:) = out_tmp;
    LEF{a} = LEF_tmp;
%     if opts.show
%         fprintf('\b|\n');
%     end
end

% ------------- ONLY PLOT GS with P<=0.05 ------------
if opts.show
    for a=1:GS.nb 
        opts_tmp = opts;
        out_tmp = out(a,:);
        
        if out_tmp(4)<=0.05
            GS_tmp = GS;
            %sort ranks in GS to find LEF
            %RI_tmp = sort(RI{a},opts_tmp.sort_type); 

            %Enrichment Score calculation
            Phit = zeros(n,1); 
            Pmiss = ones(n,1) * 1/(n-GS_tmp.entrez_nb(a));
            Phit(ind_d{a}) = (abs(RI{a}).^opts_tmp.p)/sum(abs(RI{a}));
            Pmiss(ind_d{a}) = 0;
            ES = cumsum(Phit - Pmiss);    
            [~,ind] = max(abs(ES));     %find maximum ES 
            ESmax = ES(ind);
            plot_classic(ES,ESmax,data_ranks,ind_d{a},GS.ID{a},Phit,Pmiss);
        end
        
    end
end

% ---------------------------------------------

%calculate NES p-values using NES null distribution for pos and neg
%separately
NES_perm_pos = NES_perm(NES_perm>0);
NES_perm_neg = NES_perm(NES_perm<0);
ind_pos = out(:,3) >= 0;
for a=1:GS.nb
    if ind_pos(a)
        out(a,6) = max(sum(NES_perm_pos >= out(a,5))/length(NES_perm_pos),1/length(NES_perm_pos));          %NES p-value
        out(a,7) = min(out(a,6)/(sum(out(ind_pos,5) >= out(a,5))/sum(ind_pos)),1);%NES q-value
    else
        out(a,6) = max(sum(NES_perm_neg <= out(a,5))/length(NES_perm_neg),1/length(NES_perm_neg));
        out(a,7) = min(out(a,6)/(sum(out(~ind_pos,5) <= out(a,5))/sum(~ind_pos)),1);
    end
end

out(ind_pos,8) = min(1,sum(ind_pos)*out(ind_pos,6)); % NES FWER
out(~ind_pos,8) = min(1,sum(~ind_pos)*out(~ind_pos,6));

%divide results into pos and neg
res_pos = out(ind_pos,:);
LEF_pos = LEF(ind_pos);
res_neg = out(~ind_pos,:);
LEF_neg = LEF(~ind_pos);

%sort result by GS significance
switch opts.GS_sort_type
    case 'ES_pval'
        [~,ind1] = sort(res_pos(:,4));
        [~,ind2] = sort(res_neg(:,4));
    case 'NES_qval'
        [~,ind1] = sort(res_pos(:,6));
        [~,ind2] = sort(res_neg(:,6));
    case 'NES_FWER'
        [~,ind1] = sort(res_pos(:,7));
        [~,ind2] = sort(res_neg(:,7));
    case 'none'
        ind1 = 1:size(res_pos,1);
        ind2 = 1:size(res_neg,1);
    otherwise
        error('Wrong GS sorting method selected.')
end
res_pos = res_pos(ind1,:);
ind_tmp = res_pos(:,1);
res_pos = num2cell(res_pos);
res_pos(:,1) = GS.ID(ind_tmp);

res_neg = res_neg(ind2,:);
ind_tmp = res_neg(:,1);
res_neg = num2cell(res_neg);
res_neg(:,1) = GS.ID(ind_tmp);
% res_descr={'GSNAME','GSsize','ES','ES_pval','NES','NES_pval','NES_qval','NES_FWER','NumGenes_in_LEF','LEF_genes'};
res_descr={'GSID','GSsize','ES','ES_pval','NES','NES_pval','NES_qval','NES_FWER','NumGenes_in_LEF'};

% opts.save=true;
if opts.save
    warning off
    
    % save results in .mat file
    save([ '.mat'],'res_neg','res_pos','res_descr','LEF_neg','LEF_pos')
    
    if ispc
        % save results in .xls file
        xlswrite(file_name,res_descr,'pos','A1')
        if ~isempty(res_pos)
            xlswrite(file_name,res_pos,'pos','A2');
            xlswrite(file_name,LEF_pos,'pos','J2')
        end
        xlswrite(file_name,res_descr,'neg','A1')
        if ~isempty(res_neg)
            xlswrite(file_name,res_neg,'neg','A2');
            xlswrite(file_name,LEF_neg,'neg','J2')
        end
    else
        disp('Saving to .xls on that platform is not supported.')
    end
end
