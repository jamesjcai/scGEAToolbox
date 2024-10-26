function bc = bcscore(sc, gs, expr_thres)
    if nargin < 3
        expr_thres = 0.1;
    end
    % gs - gene set    
    % gs.mode - gene set of up- and down- regulated genes
    
    % --- Checks ---
    % Check if sc is a Seurat object or expression matrix.
    if isfield(sc, 'meta.data') && isfield(sc, 'assays')
        input = 'Seurat object';
        default = sc.assays.DefaultAssay;
        if ismember(default, {'RNA', 'SCT', 'Spatial'})
            disp(['Using ', default, ' assay as input.']);
            if isequal(sc.assays.(default).counts, sc.assays.(default).data)
                error('@counts and @data matrices are identical. Is @data normalized?');
            else
                expr_matrix = sc.assays.(default).data;
            end
        else
            error('Seurat default assay must be either RNA, Spatial, or SCT.');
        end
    elseif isnumeric(sc) && ismatrix(sc)
        input = 'expression matrix';
        default = {};
        warning('Using count matrix as input. Please, check that this matrix is normalized and unscaled.');
        expr_matrix = sc;
        sc = CreateSeuratObject(expr_matrix);
    else
        error('sc must be either a Seurat object or a single-cell expression matrix.');
    end
    
    % Check if gs is a geneset object.
    if ~isa(gs, 'geneset')
        error('gs must be a geneset object.');
    end
    
    % Check expr_thres.
    if length(expr_thres) ~= 1 || expr_thres < 0 || expr_thres > 1
        error('expr.thres must be a positive number between 0 and 1.');
    end
    
    % Check gene case format
    sc_gene_case = names(CaseFraction(rownames(expr_matrix)));
    gs_gene_case = names(CaseFraction(unique(unlist({gs.genelist}))));
    if ~strcmp(sc_gene_case, gs_gene_case)
        warning(['gs genes are ', gs_gene_case, ' and sc genes are ', sc_gene_case, '. Please check your ', input, ' and translate the genes if necessary.']);
    end
    
    % --- Code ---
    % Create beyondcell object.
    bc = beyondcell(expr_matrix, sc.meta.data, sc.reductions, gs.n_genes, gs.mode, expr_thres);
    
    % Check spatial object
    if isfield(sc, 'images')
        bc.SeuratInfo.images = sc.images;
    end
    
    % Convert gene names in expr_matrix and gs to lowercase.
    expr_matrix = lower(rownames(expr_matrix));
    gs.genelist = cellfun(@(x) struct('up', lower(x.up), 'down', lower(x.down)), gs.genelist, 'UniformOutput', false);
    
    % Genes in expr_matrix.
    genes = rownames(expr_matrix);
    
    % Filter out signatures that do not share any genes with sc.
    in_sc = cellfun(@(x) any(ismember(genes, x)), gs.genelist);
    if all(~in_sc)
        error('No common genes between gs and sc were found. Stopping the execution.');
    elseif any(~in_sc)
        warning(['The following signatures do not share any genes with sc and will be removed: ', strjoin(gs.genelist(~in_sc), ', ')]);
        gs.genelist = gs.genelist(in_sc);
    end
    
    % Progress bar setup
    len_gs = length(gs.genelist);
    total = len_gs + (length(gs.mode) * len_gs);
    bins = ceil(total / 100);
    
    % Expression threshold calculation
    below_thres = false(len_gs, size(expr_matrix, 2));
    for i = 1:len_gs
        all_genes = unique(unlist(gs.genelist{i}));
        sub_expr_matrix = expr_matrix(ismember(genes, all_genes), :);
        n_expr_genes = sum(sub_expr_matrix > 0, 1);
        below_thres(i, :) = n_expr_genes < (length(all_genes) * expr_thres) | n_expr_genes == 0;
    end
    
    % Handle signatures with no passing cells
    nan_rows_idx = find(sum(below_thres, 2) == size(below_thres, 2));
    if ~isempty(nan_rows_idx)
        warning(['The following signatures have no cells that pass the expr.thres and will be removed: ', strjoin(gs.genelist(nan_rows_idx), ', ')]);
        below_thres(nan_rows_idx, :) = [];
        gs.genelist(nan_rows_idx) = [];
    end
    
    % Check remaining signatures
    if isempty(gs.genelist)
        error('No signature left. Stopping the execution.');
    end
    
    % Calculate BCS
    bcs = cell(1, length(gs.mode));
    for j = 1:length(gs.mode)
        score = zeros(len_gs, size(expr_matrix, 2));
        for k = 1:len_gs
            common_genes = intersect(genes, gs.genelist{k}.(gs.mode{j}));
            sub_expr_matrix = expr_matrix(ismember(genes, common_genes), :);
            raw = mean(sub_expr_matrix, 1);
            sig_stdev = std(sub_expr_matrix, 1, 2);
            sum_expr = sum(sub_expr_matrix, 1);
            score(k, :) = raw .* ((sum_expr - sig_stdev) ./ (raw + sig_stdev));
        end
        bcs{j} = score;
    end
    
    % Operate up and down scores
    if length(gs.mode) == 2
        nan_cells = isnan(bcs{1}) & isnan(bcs{2});
        bcs{1}(isnan(bcs{1})) = 0;
        bcs{2}(isnan(bcs{2})) = 0;
        scoring_matrix = bcs{1} - bcs{2};
        scoring_matrix(nan_cells | below_thres) = NaN;
    else
        scoring_matrix = bcs{1};
        if strcmp(gs.mode{1}, 'down')
            scoring_matrix = -scoring_matrix;
        end
        scoring_matrix(below_thres) = NaN;
    end
    
    % Invert score if necessary
    if gs.inverse_score
        scoring_matrix = -scoring_matrix;
    end
    
    bc.normalized = scoring_matrix;
    
    % Scale matrix for visualization
    bc.scaled = rescale(bc.normalized, 0, 1);
    
    % Compute switch point
    bc.switch_point = SwitchPoint(bc);
    
    disp(['There are ', num2str(sum(~isnan(bc.normalized), 'all')), '/', num2str(size(bc.normalized, 2)), ' cells without missing values in your beyondcell object.']);
end
