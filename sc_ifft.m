function [Y] = sc_ifft(X, n, k)
% https://github.com/Sanofi-Public/PMCB-scGFT/blob/master/R/Core.R
% n = nsynth            Number of synthetic cells to be generated.
% k = ncpmnts           The number of components to modify during the
%                             synthesis process.

X = pkg.norm_libsize(X, 1e4);
X = log1p(X);
ft_mtx = fft(X, [], 2);

Y = ifft(ft_mtx, [], 2);

end

%  Performs Discrete and INverse Fourier Transforms
function syn_mtx = PerformDIFT(var_mtx, groups, nsynth, ncpmnts, adj_mtx, random_seed)
    if nargin < 4
        ncpmnts = 1;
    end
    if nargin < 6
        random_seed = false;
    end
    
    if random_seed
        rng('shuffle');
    end
    
    genes_nm = rownames(var_mtx);
    cells_bc = colnames(var_mtx);
    
    % Step 1: Calculate the number of cells to synthesize per group
    groups_cnt = groupcounts(groups);
    nper_grp = round(nsynth * groups_cnt / sum(groups_cnt));
    
    excess_cells = sum(nper_grp) - nsynth;
    if excess_cells > 0
        [~, max_group] = max(nper_grp);
        nper_grp(max_group) = nper_grp(max_group) - excess_cells;
    end
    
    missing_cells = nsynth - sum(nper_grp);
    if missing_cells > 0
        [~, min_group] = min(nper_grp);
        nper_grp(min_group) = nper_grp(min_group) + missing_cells;
    end
    
    % Step 1: Apply Fourier Transform to each cell
    fprintf('Discrete Fourier transform...\n');
    ft_mtx = fft(var_mtx, [], 2);
    len_ft = size(ft_mtx, 1);          % number of genes
    
    sorted_grp_ids = sortrows(groups_cnt, 'descend');
    
    if size(ft_mtx, 2) > 1
        distr_ls = modesVar(sorted_grp_ids, cells_bc, groups, ft_mtx);
    end
    
    fprintf('Inverse Fourier transform...\n');
    fprintf('Synthesizing %d cells...\n', sum(nper_grp));
    
    synt_ls = {};
    devs_c = [];
    cell_cnt = 0;
    
    for i = 1:length(sorted_grp_ids)
        grp_id = sorted_grp_ids{i};
        cells_in_grp = cells_bc(groups == grp_id);
        
        if length(cells_in_grp) > 1
            ModificationFunc = @(cmpnts) ModificationGroup(cmpnts, distr_ls, grp_id);
        else
            ModificationFunc = @(cmpnts) ModificationSingle(cmpnts);

            % randn(length(cmpnts), 1) * 1 + 2

        end
        
        n_synt = nper_grp(i);
        cells_smpld = datasample(cells_in_grp, n_synt, 'Replace', true);
        
        synt_res = cell(1, n_synt);
        for j = 1:n_synt
            x = cells_smpld(j);
            modifiedFT = ft_mtx(:, x);
            
            cmpnts = randperm(len_ft-1, ncpmnts) + 1;
            modFactors = ModificationFunc(cmpnts);
            
            modifiedFT(cmpnts) = modifiedFT(cmpnts) .* modFactors;
            boolid = cmpnts ~= 1 & (mod(len_ft, 2) == 0 & cmpnts ~= len_ft/2 + 1);
            conjugateComps = len_ft - cmpnts(boolid) + 2;
            modifiedFT(conjugateComps) = modifiedFT(conjugateComps) .* modFactors(boolid);
            
            syn_cell = real(ifft(modifiedFT, [], 1) / len_ft);
            devs = devs_f(var_mtx(:, x), max(0, syn_cell));
            
            neighbors = find(adj_mtx(x, cells_in_grp) > 0);
            syn_cell = mean([syn_cell, var_mtx(:, neighbors)], 2);
            
            synt_res{j} = struct('synt', syn_cell, 'devs', devs);
        end
        
        synt = cellfun(@(x) x.synt, synt_res, 'UniformOutput', false);
        synt_ls = [synt_ls, synt];
        
        devs = cellfun(@(x) x.devs, synt_res);
        devs_c = [devs_c, devs];
        
        cell_cnt = cell_cnt + length(synt);
        fprintf('%d cells synthesized...\n', cell_cnt);
    end
    
    assert(length(synt_ls) == nsynth);
    syn_mtx = convert_list_to_matrix(synt_ls);
    syn_mtx(syn_mtx < 1e-6) = 0;
    
    fprintf('Deviation from originals: %.2f +/- %.2f\n', mean(devs_c), std(devs_c));
end
