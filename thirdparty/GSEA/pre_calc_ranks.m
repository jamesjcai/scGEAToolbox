function ranks_perm = pre_calc_ranks(data_ranks,data,group,opts)
% Calculate ranks for permutation test.

[m,n] = size(data);
ranks_perm = zeros(m,opts.perm_nb);
rank_type = opts.rank_type;
if_abs = opts.abs;
if_trank = opts.tied_rank;
if strcmp(opts.perm_type,'entrez')
    parfor a=1:opts.perm_nb
        ind = RandPermFast(m);
        ranks_perm(:,a) = data_ranks(ind);
    end
elseif strcmp(opts.perm_type,'pheno')
    parfor a=1:opts.perm_nb
        ind = RandPermFast(n);
        ranks_perm(:,a) = calc_ranks(data,group(ind),rank_type,if_abs,if_trank);
    end
end