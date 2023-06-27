function [opts, new_fea] = learn_LSC(fea, params)
    % learn parameters for LSC
    % input:
    %   fea: gene expression matrix
    %   params.n_clusters: number of clusters 
    %   params.HV: (yes, no) select HV
    %   params.PCA: (yes, no) apply PCA

    n_select_genes = 2000; % number of selected genes in HV.
    [m, n] = size(fea);
    n_clusters = params.n_clusters;
    % select top params.HV genes
    if (n > n_select_genes & params.HV ~= 0)
        if (params.print ~= 0)
            fprintf("Select HV genes \n");
        end
        [dump, varidx] = sort(var(fea), 'descend');
        clear dump;
        topgenes = varidx(1:n_select_genes);
        fea = fea(:,topgenes);		
    end

    % apply PCA 
    if (n > 900 & params.PCA ~= 0)
        if (params.print ~= 0)
            fprintf("Apply PCA\n");
        end
        pcaDim = 5*n_clusters;
        fea = fpca(fea, pcaDim);
        % [U, fea] = pca(fea);
        % fea = fea(:,1:pcaDim);
        clear U;
        % Write to file
        % csvwrite(strcat('/data/hoan/spectral_clustering_matlab/data/Zheng_68kPBMC-prepare-log_count_pca.csv'), [params.labels fea]);
    end

    % Input parameters
    % fprintf('use random for selection of representatives\n');
    opts.mode = 'kmeans';
    opts.r = n_clusters;
    opts.kmMaxIter = 1000;%TODO 1000
    if n < 10000 
        % opts.p= max(min(round(m/2), 20*n_clusters), round(m/5)); %round(m/10);
        opts.p= min(20*n_clusters, round(m/5)); %round(m/10);
    else 
        opts.p= min(40*n_clusters, round(m/5)); %round(m/10);
    end
    opts.maxIter = 1000;%TODO 3000
    opts.numRep = 5;%TODO 10
    opts.reduceRatio = 300;
    if isfield(params,'mode') & params.mode == 0 & m > 6000
        % fprintf('set for single mode\n');
        opts.maxIter = 1000;%
        opts.numRep = 3;%
        opts.r = max(4, n_clusters -3);
        opts.p= min(300, round(m/5)); %round(m/10);
    end

    if isfield(params,'mode') & params.mode == 2
        fprintf('set for single mode\n');
        opts.maxIter = 3000;
        opts.numRep = 10;
        opts.r = max(4, n_clusters);
        opts.p= min(500, round(m/5)); %round(m/10);
    end
    
    opts.seed = 0; % set seed
     
    % use Matlab kmeans or liteKmeans
    % opts.kmeans = 'lite';
    % factor = m/opts.reduceRatio;
    % new_fea = factor*fea/max(max(fea)); % very good for LSC Zeisel
    new_fea = fea;
end
