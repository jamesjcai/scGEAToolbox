function result = sc_nnmf(expr, rank, sel_clusters, clusterStashName, savePath)
    arguments
        expr struct          % Seurat-like object structure with 'RNA' field
        rank double = 50      % Default NMF rank
        sel_clusters double = [] % Clusters to analyze
        clusterStashName char = "default"  % Cluster stash name
        savePath char = []    % Directory to save results
    end

% Expression program signatures analysis
% We applied non-negative matrix factorization (NMF) to identify potential expression program signatures in unsupervised ways.
%In runScAnnotation, following arguments can determine detailed setting of this step.
%bool.runExprProgram indicates whether to run NMF to identify expression programs.
%nmf.rank indicates the decomposition rank used in NMF.

    % Extract gene expression data (assumed to be stored as a sparse matrix)
    data = expr.RNA.data;

    % Subset data for selected clusters if provided
    if ~isempty(sel_clusters)
        cluster_data = expr.metaData.(clusterStashName);
        data = data(:, ismember(cluster_data, sel_clusters));
    end

    % Adjust rank if larger than data dimensions
    [num_genes, num_cells] = size(data);
    if rank > min([num_genes, num_cells])
        rank = min([num_genes, num_cells]);
        warning("Input rank larger than data size, using adjusted rank.");
    end

    % Compute average gene expression
    ave_data = sum(data, 2) ./ sum(data > 0, 2);
    data = bsxfun(@minus, data, ave_data); 
    data(data < 0) = 0;

    % Perform NMF
    [W, H] = nnmf(data, rank, 'Algorithm', 'mult', 'Replicates', 10);

    % Label W and H matrices
    W_labels = arrayfun(@(x) sprintf('p%d', x), 1:size(W, 2), 'UniformOutput', false);
    H_labels = arrayfun(@(x) sprintf('p%d', x), 1:size(H, 1), 'UniformOutput', false);
    
    % Process gene weights to select top genes per program
    program_gene_value = table();
    threshold = quantile(W(:), 1 - 50 / num_genes);
    for i = 1:size(W, 2)
        gene_mask = W(:, i) > threshold;
        if sum(gene_mask) > 10
            selected_genes = W(gene_mask, i);
        else
            [~, top_genes] = maxk(W(:, i), 10);
            selected_genes = W(top_genes, i);
        end
        tmp_table = table(W_labels(i), find(gene_mask), selected_genes, 'VariableNames', {'Program', 'Gene', 'Value'});
        program_gene_value = [program_gene_value; tmp_table];
    end

    % Save results if savePath provided
    if ~isempty(savePath)
        if ~isfolder(savePath)
            mkdir(savePath);
        end
        writetable(array2table(W, 'VariableNames', W_labels), fullfile(savePath, 'W_matrix.txt'), 'Delimiter', '\t');
        writetable(array2table(H, 'RowNames', H_labels), fullfile(savePath, 'H_matrix.txt'), 'Delimiter', '\t');
        writetable(program_gene_value, fullfile(savePath, 'program_gene_values.txt'), 'Delimiter', '\t');
    end

    % Return result structure
    result.W = W;
    result.H = H;
    result.program_gene_value = program_gene_value;
end

%{
#' runExprProgram
#'
#' Perform non-negative matrix factorization (NMF) to identify expression programs.
#'
#' @param expr A Seurat object.
#' @param rank An integer of decomposition rank used in NMF.
#' @param sel.clusters A vector of selected clusters to analyze. The default is NULL and all clusters will be used.
#' @inheritParams runScAnnotation
#'
#' @return A list of decomposed matrixes (W and H), and the relative genes of each programs.
#' @export
#'
#' @importFrom methods as
#'
runExprProgram <- function(expr, rank = 50, sel.clusters = NULL, clusterStashName = "default", savePath = NULL){
    message("[", Sys.time(), "] -----: expression programs analysis")

    data <- as(object = expr[["RNA"]]@data, Class = "TsparseMatrix")
    if(!is.null(sel.clusters)){
        data <- data[, expr@meta.data[[clusterStashName]] %in% sel.clusters]
    }

    if(rank > min(dim(data))){
        rank <- min(rank, min(dim(data)))
        cat("- Warning in 'runExprProgram':
            The input rank is larger than the size of data for NMF, and use the minimum of them instead.\n")
    }

    ave.data <-  Matrix::rowSums(data) / Matrix::rowSums(data > 0)

    data@x <- data@x - ave.data[data@i + 1]
    data@x[data@x < 0] <- 0

    nmf.results <- nnmf(as.matrix(data), k = rank, verbose = 0)

    W <- nmf.results$W
    colnames(W) <- paste0("p", 1:dim(W)[2])
    H <- nmf.results$H
    rownames(H) <- paste0("p", 1:dim(H)[1])

    all.genes <- rownames(W)
    sel.W <- (W > quantile(W, 1 - 50/dim(W)[1]))
    for(pi in 1:dim(W)[2]){
        if(sum(sel.W[, pi]) > 10){
            tmp <- data.frame(program = colnames(W)[pi], gene = all.genes[sel.W[, pi]], value = W[sel.W[, pi], pi])
            tmp <- tmp[order(tmp$value, decreasing = T), ]
        }else{
            tmp <- data.frame(program = colnames(W)[pi], gene = all.genes, value = W[, pi])
            tmp <- tmp[order(tmp$value, decreasing = T), ]
            tmp <- tmp[1:10, ]
        }
        if(pi == 1){
            program.gene.value <- tmp
        }else{
            program.gene.value <- rbind(program.gene.value, tmp)
        }
    }

    # programs.geneList <- apply(sel.W, 2, FUN = function(x){ return(all.genes[x])})

    if(!is.null(savePath)){
        if(!dir.exists(file.path(savePath, "expr.programs/"))){
            dir.create(file.path(savePath, "expr.programs/"), recursive = T)
        }
        write.table(W, file = file.path(savePath, "expr.programs/W-gene-program.txt"),
                    quote = F, sep = "\t")
        write.table(H, file = file.path(savePath, "expr.programs/H-program-cell.txt"),
                    quote = F, sep = "\t")
        write.table(program.gene.value, file = file.path(savePath, "expr.programs/program.gene.value.txt"),
                    quote = F, sep = "\t", row.names = F)

        # cat("", file = file.path(savePath, "expr.programs/programs.geneList.txt"))
        # for(p in names(programs.geneList)){
        #     cat(p, "\t", str_c(programs.geneList[[p]], collapse = ", "), "\n", append = T,
        #         file = file.path(savePath, "expr.programs/programs.geneList.txt"))
        # }
    }
    return(list(W = W, H = H,
                program.gene.value = program.gene.value))
}





#' plotExprProgram
#'
#' @param H The decomposed right matrix H.
#' @param cell.annotation A data.frame of cells' annotation containing cluster information.
#' @param bool.limit A logical value indicating whether to set upper and lower limit when plot heatmap.
#' @param sel.clusters A vector of selected clusters to analyze. The default is NULL and all clusters will be used.
#' @inheritParams runScAnnotation
#'
#' @return A heatmap for cells' expression programs.
#' @export
#' @importFrom NNLM nnmf
#'
plotExprProgram <- function(H, cell.annotation, bool.limit = T, sel.clusters = NULL, savePath = NULL){
    if(bool.limit){
        up.bound <- quantile(as.matrix(H), 0.995)
        H <- limitData(H, max = up.bound)
    }
    if(!is.null(sel.clusters)){
        cell.annotation <- subset(cell.annotation, Cluster %in% sel.clusters)
    }

    tmp.results <- getClusterInfo(cell.annotation)
    cluster.info <- tmp.results$cluster.info
    cluster.colors <- tmp.results$cluster.colors
    cluster.pos <- tmp.results$cluster.pos

    p <- pheatmap(H[, rownames(cluster.info)],
                  show_colnames = F,
                  cluster_cols = F, fontsize = 7,
                  annotation_col = cluster.info,
                  annotation_colors = cluster.colors,
                  gaps_col = cluster.pos,
                  color = colorRampPalette(colors = c("#f9fcfb","#009b45"))(100),
                  silent = T)

    if(!is.null(savePath)){
        exprProgPlot.height <- 0.5 + 0.11 * dim(H)[1]
        ggsave(filename = file.path(savePath, "figures/exprProgram-heatmap.png"),
               p, width = 10, height = exprProgPlot.height, dpi = 300)
    }

    # clusters <- unique(cell.annotation$Cluster)
    # clusters <- sort(clusters)
    #
    # def.colors <- getDefaultColors(n = length(clusters))
    # cluster.colors <- c()
    # for(i in 1:length(clusters)){
    #     # cluster.colors[as.character(clusters[i])] = def.colors[clusters[i]]
    #     cluster.colors[i] = def.colors[clusters[i]]
    # }
    # cluster.colors = list(Cluster = cluster.colors)
    # ha <- HeatmapAnnotation(df = data.frame(Cluster = cell.annotation$Cluster),
    #                         name ="Cluster", col = cluster.colors)
    #
    # p <- Heatmap(H, name = "H",
    #              col = c("#f9fcfb","#009b45"),
    #              top_annotation = ha,
    #              column_split = cell.annotation$Cluster,
    #              cluster_column_slices = F,
    #              show_column_names = F,
    #              show_heatmap_legend = F)
    # if(!is.null(savePath)){
    #     png(filename = file.path(savePath, "figures/exprProgram-heatmap.png"), width = 1300, height = 800)
    #     p
    #     dev.off()
    # }

    return(p)
}
%}