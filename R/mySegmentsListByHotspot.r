## différences avec la fonction d'Isuru :
## (+) argument includeHotspots
## (-) argument minMarkers (pas besoin pour nous)
## renvoie la liste au lieu d'un objet hotspotsegments...
## un peu modifié le verbose...
## * supprimé la liste intermédiaire VII qui ne sert vraiment à rien
##   (ça pourrait être fait dans la fonction d'Isuru aussi)
## * renvoie directement le début et la fin des segments...

mySegmentsListByHotspots <- function(bedmatrix, intensity = 10, hotspots = hotspot_hg19, includeHotspots = TRUE, verbose = TRUE) {
  if(class(bedmatrix)[1] != "bed.matrix") 
    stop("Need a bed.matrix")

  if(verbose) 
    cat(paste("Using hotspots from ", deparse(substitute(hotspots)), "\n"))

  dataFrameColNames <- c("Chromosome", "Start", "End", "IntensitycMMb")
  
  if (!all(dataFrameColNames %in% colnames(hotspots))) 
    stop("'hotspots' should have columns names Chromosome, Start, End, IntensitycMMb")
  
  if(verbose) 
    cat("Listing hotspots\n")
  
  chr.ids <- as.character(intersect(unique(bedmatrix@snps$chr), unique(hotspots$Chromosome)))
  VI <- list()
 
  if(verbose) cat("Chr")
  for (i in chr.ids) {
    if (verbose) cat(" ",i)
    chr_hotspot <- hotspots[which(hotspots$Chromosome == i), ]
    w <- which(chr_hotspot$IntensitycMMb > intensity)
    if(includeHotspots)
      segment <- cbind(c(0, chr_hotspot$Centre[w])+1, c(chr_hotspot$Centre[w], Inf))
    else
      segment <- cbind(c(0, chr_hotspot$End[w]), c(chr_hotspot$Start[w], Inf))
    VI[[i]] <- segment
  }
  if (verbose) cat("\n")
  nb.hs <- sum(sapply(VI,nrow))
  if (verbose) cat("Found ", nb.hs, " segments between hotspots with intensity > ", intensity, "\n")
  
  if (verbose) cat("Distributing markers in segments...\n")

  shift <- sapply(chr.ids, function(i) which(bedmatrix@snps$chr == i)[1]) - 1L

  debut <- integer(nb.hs)
  fin   <- debut
  k <- 1
  for (i in chr.ids) {
    chr_segment <- VI[[i]]
    mkr <- bedmatrix@snps$pos[bedmatrix@snps$chr == i]
    chr <- list()
    for (j in 1:nrow(chr_segment)) {
      b <- which(mkr > chr_segment[j, 1] & mkr < chr_segment[j, 2])
      if (length(b) == 0) next
      debut[k] <- shift[i] + b[1]
      fin[k]   <- shift[i] + b[length(b)]
      k <- k+1
    }
  }
  debut <- debut[1:(k-1)]
  fin <- fin[1:(k-1)]
  list(debut = debut, fin = fin)
}

