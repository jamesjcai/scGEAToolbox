function s = sc_stemsig(X, stem_sig, species)
    arguments
        X double {mustBeNonempty}  % Ensure the gene expression matrix is a non-empty numeric array
        stem_sig double = []       % Default is empty, optional argument
        species (1,:) char {mustBeMember(species, ["human", "mouse"])} = "human"  % Either 'human' or 'mouse'
    end    
    
    % Load default stemness signature if not provided
    if isempty(stem_sig)
        stem_sig = readtable('pcbc-stemsig.tsv', 'FileType', 'text', 'ReadVariableNames', false, 'RowNames', true);
        if strcmp(species, 'mouse')
            % Convert human genes to mouse equivalents (assuming getMouseGene is implemented)
            sig_genes = getMouseGene(stem_sig.Properties.RowNames, true);
            stem_sig = stem_sig(sig_genes, :);
            stem_sig.Properties.RowNames = sig_genes;
        end
    end
    
    % Ensure common genes between stem signature and input matrix
    common_genes = intersect(stem_sig.Properties.RowNames, X.Properties.RowNames);
    X_common = X(common_genes, :);
    stem_sig_common = stem_sig(common_genes, :);
    
    % Initialize the stemness score array
    s = zeros(1, width(X_common));
    
    % Calculate Spearman correlation for each cell (column in X)
    for i = 1:width(X_common)
        s(i) = corr(X_common{:, i}, stem_sig_common{:, 1}, 'Type', 'Spearman', 'Rows', 'complete');
    end
    
    % Normalize the stemness scores between 0 and 1
    s = s - min(s);
    s = s / max(s);
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