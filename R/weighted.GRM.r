#' Weighted GRM
#'
#' @param x a bed matrix
#' @param weights a vector of SNP weights
#' @param which.snps a logical vector 
#' @param autosome.only TRUE to keep only autosomes
#'
#' @return A weighted GRM
#' 
#' @export
#'
#' @examples 
#' # see example in SNP.weights
weighted.GRM <- function(x, weights = rep(1,ncol(x)), which.snps, autosome.only = TRUE) {
  if(missing(which.snps)) which.snps <- rep(TRUE, ncol(x))
  if(autosome.only) 
    which.snps <- which.snps & is.autosome(x@snps$chr)

  which.snps <- which.snps & (x@p > 0) & (x@p < 1) & (weights != 0)
  if(!is.logical(which.snps) | length(which.snps) != ncol(x))
    stop("which.snps must be a Logical vector of length ncol(x)")
  
  weights <- weights/sum(weights)
  w <- weights/(2*x@p*(1-x@p)) 
  K <- .Call('_migou_weighted_Kinship_w', x@bed, 2*x@p[which.snps], w[which.snps], which.snps, 1L) 

  if(!is.null(x@ped$id)) {
    if(anyDuplicated(x@ped$id) == 0) 
      rownames(K) <- colnames(K) <- x@ped$id
  }

  K
}
