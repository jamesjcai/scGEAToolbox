% 
% % https://rdrr.io/bioc/scran/src/R/mnnCorrect.R
% % https://github.com/yycunc/SMNN/blob/master/R/SMNNcorrect.R
% 
% x = 1:100;
% A = cos(2*pi*0.05*x+2*pi*rand) + 0.5*randn(1,100);
% [B, window] = smoothdata(A,'gaussian');
% window
% 
% Expected type of input data
% The input expression values should generally be log-transformed, e.g., log-counts, see normalize
% for details. They should also be normalized within each data set to remove cell-specific biases in
% capture efficiency and sequencing depth. By default, a further cosine normalization step is performed on the supplied expression data to eliminate gross scaling differences between data sets.
% • When cos.norm.in=TRUE, cosine normalization is performed on the matrix of expression
% values used to compute distances between cells. This can be turned off when there are no
% scaling differences between data sets.
% • When cos.norm.out=TRUE, cosine normalization is performed on the matrix of values used
% to calculate correction vectors (and on which those vectors are applied). This can be turned
% off to obtain corrected values on the log-scale, similar to the input data.
% The cosine normalization is achieved using the cosineNorm function.
% Users should note that the order in which batches are corrected will affect the final results. The
% first batch in order is used as the reference batch against which the second batch is corrected.
% Corrected values of the second batch are added to the reference batch, against which the third batch
% is corrected, and so on. This strategy maximizes the chance of detecting sufficient MNN pairs for
% stable calculation of correction vectors in subsequent batches.
% Further options
% The function depends on a shared biological manifold, i.e., one or more cell types/states being
% present in multiple batches. If this is not true, MNNs may be incorrectly identified, resulting in
% over-correction and removal of interesting biology. Some protection can be provided by removing
% components of the correction vectors that are parallel to the biological subspaces in each batch. The
% biological subspace in each batch is identified with a SVD on the expression matrix, using either
% svd or irlba. The number of dimensions of this subspace can be controlled with svd.dim. (By
% default, this option is turned off by setting svd.dim=0.)
% If var.adj=TRUE, the function will adjust the correction vector to equalize the variances of the
% two data sets along the batch effect vector. In particular, it avoids “kissing” effects whereby MNN
% pairs are identified between the surfaces of point clouds from different batches. Naive correction
% would then bring only the surfaces into contact, rather than fully merging the clouds together. The
% adjustment ensures that the cells from the two batches are properly intermingled after correction.
% This is done by identifying each cell’s position on the correction vector, identifying corresponding
% quantiles between batches, and scaling the correction vector to ensure that the quantiles are matched
% after correction.
%     
%     
%     
%     mnnCorrect <- function(..., k=20, sigma=0.1, cos.norm.in=TRUE, cos.norm.out=TRUE, 
%     
%                        svd.dim=0L, var.adj=TRUE, compute.angle=FALSE,
%                        subset.row=NULL, order=NULL, pc.approx=FALSE, irlba.args=list(),
%                        BNPARAM=KmknnParam(), BPPARAM=SerialParam())
% # Performs batch correction on multiple matrices of expression data,
% # as specified in the ellipsis.
% #    
% # written by Laleh Haghverdi
% # with modifications by Aaron Lun
% # created 7 April 2017
% {
%     .Deprecated("batchelor::mnnCorrect")
%     batches <- list(...) 
%     nbatches <- length(batches) 
%     if (nbatches < 2L) { 
%         stop("at least two batches must be specified") 
%     }
%     
%     prep.out <- prepare.input.data(batches, cos.norm.in=cos.norm.in, cos.norm.out=cos.norm.out, subset.row=subset.row)
%     in.batches <- prep.out$In
%     out.batches <- prep.out$Out
%     subset.row <- prep.out$Subset
%     same.set <- prep.out$Same
% 
%     # Setting up the order.
%     if (is.null(order)) {
%         order <- seq_len(nbatches)
%     } else {
%         order <- as.integer(order)
%         if (!identical(seq_len(nbatches), sort(order))) { 
%             stop(sprintf("'order' should contain values in 1:%i", nbatches))
%         }
%     }
%    
%     # Setting up the variables.
%     ref <- order[1]
%     ref.batch.in <- t(in.batches[[ref]])
%     if (!same.set) { 
%         ref.batch.out <- t(out.batches[[ref]])
%     }
%     output <- vector("list", nbatches)
%     output[[ref]] <- out.batches[[ref]]
%     mnn.list <- vector("list", nbatches)
%     mnn.list[[ref]] <- DataFrame(current.cell=integer(0), other.cell=integer(0), other.batch=integer(0))
%     original.batch <- rep(ref, nrow(ref.batch.in)) 
%     angle.list <- vector("list", nbatches)
%     angle.list[[ref]] <- numeric(0)
% 
%     # Looping through the batches.
%     for (b in 2:nbatches) { 
%         target <- order[b]
%         other.batch.in.untrans <- in.batches[[target]]
%         other.batch.in <- t(other.batch.in.untrans)
%         if (!same.set) { 
%             other.batch.out.untrans <- out.batches[[target]]
%             other.batch.out <- t(other.batch.out.untrans)
%         }
%         
%         # Finding pairs of mutual nearest neighbours.
%         sets <- find.mutual.nn(ref.batch.in, other.batch.in, k1=k, k2=k, BNPARAM=BNPARAM, BPPARAM=BPPARAM)
%         s1 <- sets$first
%         s2 <- sets$second      
% 
%         # Computing the correction vector.
%         correction.in <- compute.correction.vectors(ref.batch.in, other.batch.in, s1, s2, other.batch.in.untrans, sigma)
%         if (!same.set) {
%             correction.out <- compute.correction.vectors(ref.batch.out, other.batch.out, s1, s2, other.batch.in.untrans, sigma)
%             # NOTE: use of 'other.batch.in.untrans' here is intentional, 
%             # as the distances should be the same as the MNN distances.
%         }
% 
%         # Calculating the smallest angle between each correction vector and the first 2 basis vectors of the reference batch.
%         if (compute.angle) {
%             ref.centred <- t(ref.batch.in)
%             ref.centred <- ref.centred - rowMeans2(DelayedArray(ref.centred))
%             ref.basis <- .svd_internal(ref.centred, nu=2, nv=0, pc.approx=pc.approx, irlba.args=irlba.args)$u
% 
%             angle.out <- numeric(nrow(correction.in))
%             for (i in seq_along(angle.out)) {
%                 angle.out[i] <- find.shared.subspace(ref.basis, t(correction.in[i,,drop=FALSE]))$angle
%             }
%             angle.list[[target]] <- angle.out
%         }
% 
%         # Removing any component of the correction vector that's parallel to the biological basis vectors in either batch.
%         # Only using cells involved in MNN pairs, to avoid undercorrection in directions that weren't problematic anyway
%         # (given that using SVDs was intended to mitigate the effect of identifying the wrong MNN pairs).
%         if (svd.dim>0) {
%             u1 <- unique(s1)
%             u2 <- unique(s2)
% 
%             # Computing the biological subspace in both batches, and subtract it from the batch correction vector.
%             in.span1 <- get.bio.span(t(ref.batch.in[u1,,drop=FALSE]), ndim=svd.dim,
%                                      pc.approx=pc.approx, irlba.args=irlba.args)
%             in.span2 <- get.bio.span(other.batch.in.untrans[,u2,drop=FALSE], ndim=svd.dim, 
%                                      pc.approx=pc.approx, irlba.args=irlba.args)
%             correction.in <- subtract.bio(correction.in, in.span1, in.span2)
% 
%             # Repeating for the output values.
%             if (!same.set) { 
%                 out.span1 <- get.bio.span(t(ref.batch.out[u1,,drop=FALSE]), subset.row=subset.row,
%                                           ndim=svd.dim, pc.approx=pc.approx, irlba.args=irlba.args)
%                 out.span2 <- get.bio.span(other.batch.out.untrans[,u2,drop=FALSE], subset.row=subset.row,
%                                           ndim=svd.dim, pc.approx=pc.approx, irlba.args=irlba.args)
%                 correction.out <- subtract.bio(correction.out, out.span1, out.span2, subset.row=subset.row)
%             }
%         } 
%        
%         # Adjusting the shift variance; done after any SVD so that variance along the correction vector is purely technical.
%         if (var.adj) { 
%             correction.in <- adjust.shift.variance(ref.batch.in, other.batch.in, correction.in, sigma=sigma)
%             if (!same.set) {
%                 correction.out <- adjust.shift.variance(ref.batch.out, other.batch.out, correction.out, sigma=sigma, subset.row=subset.row) 
%             }
%         }
% 
%         # Applying the correction and expanding the reference batch. 
%         other.batch.in <- other.batch.in + correction.in
%         ref.batch.in <- rbind(ref.batch.in, other.batch.in)
%         if (same.set) {
%             output[[target]] <- t(other.batch.in)
%         } else {
%             other.batch.out <- other.batch.out + correction.out
%             ref.batch.out <- rbind(ref.batch.out, other.batch.out)
%             output[[target]] <- t(other.batch.out)
%         }
% 
%         # Storing the identities of the MNN pairs (using RLEs for compression of runs).
%         mnn.list[[target]] <- DataFrame(current.cell=s2, other.cell=Rle(s1), other.batch=Rle(original.batch[s1]))
%         original.batch <- c(original.batch, rep(target, nrow(other.batch.in)))
%     }
% 
%     # Formatting output to be consistent with input.
%     names(output) <- names(batches)
%     names(mnn.list) <- names(batches)
%     final <- list(corrected=output, pairs=mnn.list)
%     if (compute.angle){ 
%         names(angle.list) <- names(batches)
%         final$angles <- angle.list
%     }
%     return(final)
% }
% 
% 
