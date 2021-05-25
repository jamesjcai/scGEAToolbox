function ES_perm = ES_perm_test(ranks_perm,data_entrez,GS_entrez,perm_nb,sort_type)
% Permutation test to find distribution of Enrichment Score

n = length(data_entrez);
ES_perm = zeros(1,perm_nb);
GS_entrez_nb = length(GS_entrez);
for a=1:perm_nb
    [data_ranks,idx] = sort(ranks_perm(:,a),sort_type);
    data_entrez_tmp = data_entrez(idx);

    %find genes with ranks within GS
    [~,ind_d] = intersect_fast(data_entrez_tmp,GS_entrez);  %analyzed genes that are in given GS

    %Enrichment Score calculation
    Phit = zeros(n,1); 
    Pmiss = ones(n,1) * 1/(n-GS_entrez_nb);
    Phit(ind_d) = abs(data_ranks(ind_d))/sum(abs(data_ranks(ind_d)));
    Pmiss(ind_d) = 0;
    ES = cumsum(Phit - Pmiss);
    [~,ind] = max(abs(ES));    %find maximum ES 
    ES_perm(a) = ES(ind);
end