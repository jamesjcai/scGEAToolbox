% https://satijalab.org/seurat/v3.0/sctransform_vignette.html
% 
% #' Use regularized negative binomial regression to normalize UMI count data
% #'
% #' This function calls sctransform::vst. The sctransform package is available at
% #' https://github.com/ChristophH/sctransform.
% #' Use this function as an alternative to the NormalizeData,
% #' FindVariableFeatures, ScaleData workflow. Results are saved in a new assay
% #' (named SCT by default) with counts being (corrected) counts, data being log1p(counts),
% #' scale.data being pearson residuals; sctransform::vst intermediate results are saved
% #' in misc slot of new assay.
% #'
% #' @param object A seurat object
% #' @param assay Name of assay to pull the count data from; default is 'RNA'
% #' @param new.assay.name Name for the new assay containing the normalized data
% #' @param do.correct.umi Place corrected UMI matrix in assay counts slot; default is TRUE
% #' @param ncells Number of subsampling cells used to build NB regression; default is NULL
% #' @param variable.features.n Use this many features as variable features after
% #' ranking by residual variance; default is 3000
% #' @param variable.features.rv.th Instead of setting a fixed number of variable features,
% #' use this residual variance cutoff; this is only used when \code{variable.features.n}
% #' is set to NULL; default is 1.3
% #' @param vars.to.regress Variables to regress out in a second non-regularized linear
% #' regression. For example, percent.mito. Default is NULL
% #' @param do.scale Whether to scale residuals to have unit variance; default is FALSE
% #' @param do.center Whether to center residuals to have mean zero; default is TRUE
% #' @param clip.range Range to clip the residuals to; default is \code{c(-sqrt(n/30), sqrt(n/30))},
% #' where n is the number of cells
% #' @param conserve.memory If set to TRUE the residual matrix for all genes is never
% #' created in full; useful for large data sets, but will take longer to run;
% #' this will also set return.only.var.genes to TRUE; default is FALSE
% #' @param return.only.var.genes If set to TRUE the scale.data matrices in output assay are
% #' subset to contain only the variable genes; default is TRUE
% #' @param seed.use Set a random seed. By default, sets the seed to 1448145. Setting
% #' NULL will not set a seed.
% #' @param verbose Whether to print messages and progress bars
% #' @param ... Additional parameters passed to \code{sctransform::vst}
% #'
% #' @return Returns a Seurat object with a new assay (named SCT by default) with
% #' counts being (corrected) counts, data being log1p(counts), scale.data being
% #' pearson residuals; sctransform::vst intermediate results are saved in misc
% #' slot of the new assay.
% #'
% #' @importFrom stats setNames
% #' @importFrom sctransform vst get_residual_var get_residuals correct_counts
% #'
% #' @export
% #'
% #' @examples
% #' SCTransform(object = pbmc_small)
% #'
% SCTransform <- function(
%   object,
%   assay = 'RNA',
%   new.assay.name = 'SCT',
%   do.correct.umi = TRUE,
%   ncells = NULL,
%   variable.features.n = 3000,
%   variable.features.rv.th = 1.3,
%   vars.to.regress = NULL,
%   do.scale = FALSE,
%   do.center = TRUE,
%   clip.range = c(-sqrt(x = ncol(x = object[[assay]]) / 30), sqrt(x = ncol(x = object[[assay]]) / 30)),
%   conserve.memory = FALSE,
%   return.only.var.genes = TRUE,
%   seed.use = 1448145,
%   verbose = TRUE,
%   ...
% ) {
%   if (!is.null(x = seed.use)) {
%     set.seed(seed = seed.use)
%   }
%   assay <- assay %||% DefaultAssay(object = object)
%   assay.obj <- GetAssay(object = object, assay = assay)
%   umi <- GetAssayData(object = assay.obj, slot = 'counts')
%   cell.attr <- slot(object = object, name = 'meta.data')
%   vst.args <- list(...)
%   # check for batch_var in meta data
%   if ('batch_var' %in% names(x = vst.args)) {
%     if (!(vst.args[['batch_var']] %in% colnames(x = cell.attr))) {
%       stop('batch_var not found in seurat object meta data')
%     }
%   }
%   # check for latent_var in meta data
%   if ('latent_var' %in% names(x = vst.args)) {
%     known.attr <- c('umi', 'gene', 'log_umi', 'log_gene', 'umi_per_gene', 'log_umi_per_gene')
%     if (!all(vst.args[['latent_var']] %in% c(colnames(x = cell.attr), known.attr))) {
%       stop('latent_var values are not from the set of cell attributes sctransform calculates by default and cannot be found in seurat object meta data')
%     }
%   }
%   # check for vars.to.regress in meta data
%   if (any(!vars.to.regress %in% colnames(x = cell.attr))) {
%     stop('problem with second non-regularized linear regression; not all variables found in seurat object meta data; check vars.to.regress parameter')
%   }
%   if (any(c('cell_attr', 'show_progress', 'return_cell_attr', 'return_gene_attr', 'return_corrected_umi') %in% names(x = vst.args))) {
%     warning(
%       'the following arguments will be ignored because they are set within this function:',
%       paste(
%         c(
%           'cell_attr',
%           'show_progress',
%           'return_cell_attr',
%           'return_gene_attr',
%           'return_corrected_umi'
%         ),
%         collapse = ', '
%       ),
%       call. = FALSE,
%       immediate. = TRUE
%     )
%   }
%   vst.args[['umi']] <- umi
%   vst.args[['cell_attr']] <- cell.attr
%   vst.args[['show_progress']] <- verbose
%   vst.args[['return_cell_attr']] <- TRUE
%   vst.args[['return_gene_attr']] <- TRUE
%   vst.args[['return_corrected_umi']] <- do.correct.umi
%   vst.args[['n_cells']] <- ncells
%   residual.type <- vst.args[['residual_type']] %||% 'pearson'
%   res.clip.range <- vst.args[['res_clip_range']] %||% c(-sqrt(x = ncol(x = umi)), sqrt(x = ncol(x = umi)))
%   if (conserve.memory) {
%     return.only.var.genes <- TRUE
%   }
%   if (conserve.memory) {
%     vst.args[['residual_type']] <- 'none'
%     vst.out <- do.call(what = 'vst', args = vst.args)
%     feature.variance <- get_residual_var(
%       vst_out = vst.out,
%       umi = umi,
%       residual_type = residual.type,
%       res_clip_range = res.clip.range
%     )
%     vst.out$gene_attr$residual_variance <- NA_real_
%     vst.out$gene_attr[names(x = feature.variance), 'residual_variance'] <- feature.variance
%   } else {
%     vst.out <- do.call(what = 'vst', args = vst.args)
%     feature.variance <- setNames(
%       object = vst.out$gene_attr$residual_variance,
%       nm = rownames(x = vst.out$gene_attr)
%     )
%   }
%   if (verbose) {
%     message('Determine variable features')
%   }
%   feature.variance <- sort(x = feature.variance, decreasing = TRUE)
%   if (!is.null(x = variable.features.n)) {
%     top.features <- names(x = feature.variance)[1:min(variable.features.n, length(x = feature.variance))]
%   } else {
%     top.features <- names(x = feature.variance)[feature.variance >= variable.features.rv.th]
%   }
%   if (verbose) {
%     message('Set ', length(x = top.features), ' variable features')
%   }
%   if (conserve.memory) {
%     # actually get the residuals this time
%     if (verbose) {
%       message("Return only variable features for scale.data slot of the output assay")
%     }
%     vst.out$y <- get_residuals(
%       vst_out = vst.out,
%       umi = umi[top.features, ],
%       residual_type = residual.type,
%       res_clip_range = res.clip.range
%     )
%     if (do.correct.umi & residual.type == 'pearson') {
%       vst.out$umi_corrected <- correct_counts(
%         x = vst.out,
%         umi = umi,
%         show_progress = verbose
%       )
%     }
%   }
%   # create output assay and put (corrected) umi counts in count slot
%   if (do.correct.umi & residual.type == 'pearson') {
%     if (verbose) {
%       message('Place corrected count matrix in counts slot')
%     }
%     assay.out <- CreateAssayObject(counts = vst.out$umi_corrected)
%   } else {
%     assay.out <- CreateAssayObject(counts = umi)
%   }
%   # set the variable genes
%   VariableFeatures(object = assay.out) <- top.features
%   # put log1p transformed counts in data
%   assay.out <- SetAssayData(
%     object = assay.out,
%     slot = 'data',
%     new.data = log1p(x = GetAssayData(object = assay.out, slot = 'counts'))
%   )
%   if (return.only.var.genes & !conserve.memory) {
%     scale.data <- vst.out$y[top.features, ]
%   } else {
%     scale.data <- vst.out$y
%   }
%   # clip the residuals
%   scale.data[scale.data < clip.range[1]] <- clip.range[1]
%   scale.data[scale.data > clip.range[2]] <- clip.range[2]
%   # 2nd regression
%   scale.data <- ScaleData(
%     scale.data,
%     features = NULL,
%     vars.to.regress = vars.to.regress,
%     latent.data = cell.attr[, vars.to.regress, drop = FALSE],
%     model.use = 'linear',
%     use.umi = FALSE,
%     do.scale = do.scale,
%     do.center = do.center,
%     scale.max = Inf,
%     block.size = 750,
%     min.cells.to.block = 3000,
%     verbose = verbose
%   )
%   assay.out <- SetAssayData(
%     object = assay.out,
%     slot = 'scale.data',
%     new.data = scale.data
%   )
%   # save vst output (except y) in @misc slot
%   vst.out$y <- NULL
%   Misc(object = assay.out, slot = 'vst.out') <- vst.out
%   # also put gene attributes in meta.features
%   assay.out[[paste0('sct.', names(x = vst.out$gene_attr))]] <- vst.out$gene_attr
%   assay.out[['sct.variable']] <- rownames(x = assay.out[[]]) %in% top.features
%   object[[new.assay.name]] <- assay.out
%   if (verbose) {
%     message(paste("Set default assay to", new.assay.name))
%   }
%   DefaultAssay(object = object) <- new.assay.name
%   object <- LogSeuratCommand(object = object)
%   return(object)
% }
% 
