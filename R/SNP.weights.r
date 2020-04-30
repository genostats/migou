#' Compute weights
#'
#' @param x a bed matrix
#' @param hotspots data frame of hotspots positions
#' @param hotspots.intensity min intensity of hotspot 
#' @param lambda regularisation parameter 
#' @param verbose verbosity
#' 
#' @details Computes weights on segments separated by hotspots.
#'
#' @return A SNP map with weights...
#' 
#'
#' @examples
#' # loading data (503 european individuals from 1KG project x 536k SNPs)
#' KG <- read.bed.matrix(system.file("extdata", "1KG_genos.bed", package="migou"))
#' # computing weights
#' map <- SNP.weights(KG, hotspots = hotspots.hg19)
#' # computing the weighted GRM
#' K <- weighted.GRM(KG, map$weight)
#' 
#' @export
SNP.weights <- function(x, hotspots, hotspots.intensity = 10,  lambda = 0.01, verbose = TRUE) {
    # segments par hotspots
    s <- mySegmentsListByHotspots(x, intensity = hotspots.intensity, hotspots = hotspots, verbose = verbose)
    debut <- s$debut
    fin <- s$fin
    
    nb.snps <- (fin - debut + 1)  # à comparer avec ncol(x)
    if(verbose) {
      cat(length(debut), "segments of size from", min(nb.snps), "to", max(nb.snps), "\n")
      cat("median size", median(nb.snps), "\n")
      cat("Total : ", sum(nb.snps), "SNPs\n")
    }
    
    # weights.bloc :
    # fonction qui bouffe deux coordonnées sur x, calcule la matrice de LD, 
    # pour renvoyer w + variance de l'estimateur sur le bloc
    # plus v0 = variance de l'estimateur avec poids constants
    
    ## ON LAISSE LD() standardiser x avec "mu_sigma" sinon ça
    ## part en vrille quand il y a des eigen values < 0 !!!!
    ## Ceci laisse penser qu'utiliser une estimation par maximum 
    ## de vraisemblance serait une erreur car on peut avoir le même
    ## problème
    standardize(x) <- "none"
    
    if(verbose) cat("Variance optimization on each segment\n")
    k <- 0;
    kk <- options("width")[[1]] - 6
    weights.bloc <- function(be, en) {
      if(verbose) {
        cat("."); k <<- k+1; if( (k%%kk) == 0) cat(k, "\n")
      }
      p <- x@p[be:en]
      which.non.mono <- (p > 0 & p < 1)
      
      A <- LD(x, c(be, en)) + lambda*diag(en-be+1)
      A <- A[which.non.mono, which.non.mono, drop = FALSE]
      w <- rep(1, ncol(A)) / ncol(A)
      
      # "constant weigths" variance
      v0 <- sum(w * (A %*% w))
      if(ncol(A) > 2)
        v <- min_var_eigen(w, A, ncol(A) - 1)
      else
        v <- sum(w * (A %*% w))
      
      w1 <- numeric(en-be+1)
      w1[which.non.mono] <- w
      
      list(var0 = v0, var = v, weights = w1)
    }
    # hop on l'applique
    R <- mapply(weights.bloc, debut, fin)
    cat("\n");
    
    SEG <- data.frame(chr = x@snps$chr[debut],
                      pos.debut = x@snps$pos[debut],
                      pos.fin = x@snps$pos[fin],
                      i.debut = debut,
                      i.fin = fin)
    
    # Il y a un petit paquet de SNPs non pris en compte car ils sont "dans les hotspots"
    SEG$nb.snps <- nb.snps
    
    SEG$var0 <- unlist(R["var0", ])
    SEG$var <- unlist(R["var", ])
    SEG$weights <- R["weights",]
    
    cat("Building map with weights\n")
    # On calcule les poids sur la carte de x (par ailleurs utilisée pour récupérer les haplotypes)
    MAP <- x@snps[,1:6]
    
    MAP$weight <- 0
    for(i in 1:nrow(SEG))
      MAP$weight[ debut[i]:fin[i] ] <- SEG$weights[[i]] / SEG$var[i]
    
    MAP$weight <- MAP$weight/sum(MAP$weight)
    
    # list(SEG = SEG, MAP = MAP, var =  1/sum(1/SEG$var),
    #     var0 = sum( SEG$nb.snps**2 * SEG$var0 ) / sum(SEG$nb.snps)**2)
    MAP
}

