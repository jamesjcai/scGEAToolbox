function s = sc_stemness(X, g)
    %arguments
    %    X double {mustBeNonempty}  % Ensure the gene expression matrix is a non-empty numeric array
    %end
   % https://www.nature.com/articles/sdata201730
   % https://github.com/czythu/scCancer
   % https://github.com/wguo-research/scCancer/tree/master/inst/txt
   % https://academic.oup.com/bib/article/22/3/bbaa127/5867555?login=false

pw1 = fileparts(mfilename('fullpath'));
dbfile1 = fullfile(pw1, 'resources', 'scCancer', 'pcbc_stemsig.txt');
if ~exist(dbfile1, 'file'), error('Missing file pcbc_stemsig.txt.'); end
   
    % Load default stemness signature if not provided
    T = readtable(dbfile1, 'FileType', 'text', 'ReadVariableNames', false);
        
    % Ensure common genes between stem signature and input matrix
    [~, ix, iy] = intersect(T.Var1, g);
    %X = sc_norm(X);
    %X = log1p(X);
    X = X(iy, :);
    if issparse(X), X = full(X); end
    stem_sig_common = T.Var2(ix, :);
    
    % Initialize the stemness score array
    % s = zeros(size(X, 2), 1);
    
    % Calculate Spearman correlation for each cell (column in X)
    % for i = 1:width(X_common)
    %     s(i) = corr(X_common(:, i), stem_sig_common(:), ...
    %         'Type', 'Spearman');
    % end
    
    s = corr(X, stem_sig_common, ...
        "Type","Spearman");
    
    % Normalize the stemness scores between 0 and 1
    %s = s - min(s);
    %s = s / max(s);
    s = normalize(s, "range");
end

% https://github.com/wguo-research/scCancer/blob/master/R/scAnnotation.R#L1058
%{
#' runStemness
#'
#' Estimate cell stemness according to the Spearman correlation with stemness signature.
#'
#' @param X An expression matrix of gene by cell to estimate stemness.
#' @param stem.sig An array of stemness signature. The default is NULL, and a prepared signature will be used.
#' @inheritParams runScAnnotation
#'
#' @return An array of cell stemness scores.
#' @export
#'
runStemness <- function(X, stem.sig = NULL, species = "human"){
    message("[", Sys.time(), "] -----: stemness score calculation")
    if(is.null(stem.sig)){
        stem.sig.file <- system.file("txt", "pcbc-stemsig.tsv", package = "scCancer")
        stem.sig <- read.delim(stem.sig.file, header = FALSE, row.names = 1)
        if(species == "mouse"){
            sig.genes <- getMouseGene(rownames(stem.sig), bool.name = T)
            stem.sig <- stem.sig[names(sig.genes), , drop=F]
            rownames(stem.sig) <- sig.genes
        }
    }

    common.genes <- intersect(rownames(stem.sig), rownames(X))
    X <- X[common.genes, ]
    stem.sig <- stem.sig[common.genes, ]

    s <- apply(X, 2, function(z) {cor(z, stem.sig, method = "sp", use = "complete.obs")})
    names(s) <- colnames(X)

    s <- s - min(s)
    s <- s / max(s)

    return(s)
}
%}