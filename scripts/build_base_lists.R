 
## for this script we are building the incidence matrix that we are going to use
## to improve the analysis

rm(list = ls())

library(GenomicRanges)
library(data.table)

list_dir <- "inst/base_list"

files <- list.files(list_dir)
files <- files[grep("Rename",files,invert = TRUE)]

tables <- lapply(file.path(list_dir,files),read.table,stringsAsFactors = FALSE,header = TRUE)
tables <- lapply(tables,data.table)
names(tables) <- do.call(c,lapply(strsplit(files,"_"),function(x)x[[1]]))

## remove the peak of an outlier profile
tables[["EBNA3C"]] <- tables[["EBNA3C"]][-which(chrID == "chr2" & peakStart == 33141400)]

mosaics <- lapply(tables,function(x)x[,1:6,with = FALSE])

dt2gr <- function(x)
{
  out <- GRanges(seqnames = x[,(chrID)],
                 ranges = IRanges(start = x[,(peakStart)],
                   end = x[,(peakStop)]),strand = "*")
  out$aveP <- x[,(aveP)]
  out$minP <- x[,(minP)]
  return(out)
}

mosaics_gr <- lapply(mosaics,dt2gr)

unify_list <- function(...)
{
  peak_list <- list(...)
  if(!is.list(peak_list[[1]])){
    peak_list <- list(...)
  }else{
    peak_list <- peak_list[[1]]
  }
  ov <- findOverlaps(peak_list[["JK234"]] , peak_list[["JK92"]])
  ## build RBPJ as JK234 + JK92 (that don't overlap JK234)
  ## if we want to build the peaks with JK separate and then make RBPJ as any of both
  ## all_peaks$RBPJ <- ifelse(all_peaks$JK92 + all_peaks$JK234 > 0 , 1 , 0)

  peak_list[["RBPJ"]] <- sort(c(peak_list[["JK234"]],peak_list[["JK92"]][-subjectHits(ov)]))

  
  all_peaks <- Reduce(c,peak_list[grep("JK",names(peak_list),invert = TRUE)])
  all_peaks <- reduce(all_peaks)

  overlaps <- lapply(peak_list, function(x)ifelse(countOverlaps(all_peaks,x) > 0 , 1, 0))
  overlaps <- as.data.table(overlaps)
  ov_mats <- lapply(peak_list,findOverlaps,all_peaks)

  post_prob <- mapply(function(ov_mat,gr){
    ov_mat <- data.table(as.data.frame(ov_mat))
    ov_mat[, minP := gr[ov_mat[,(queryHits)]]$minP]
    minP <- ov_mat[,min(minP),by = subjectHits][,(V1)]
    return(minP)
    },ov_mats,peak_list,SIMPLIFY = FALSE)

  probs <- copy(overlaps)
  probs <- mapply( function(x,y){
    x[x == 1] <- y
    x[x == 0] <- NA
    return(x)},probs,post_prob,SIMPLIFY  = FALSE)
  probs <- as.data.table(probs)
  probs[,minProb := pmin(EBNA2, EBNA3A, EBNA3B, EBNA3C,RBPJ,na.rm = TRUE)]

  peaks <- data.table(seqnames = as.character(seqnames(all_peaks)),start = start(all_peaks),
                end = end(all_peaks),width = width(all_peaks))

  peak_coordinates <- lapply(peak_list[grep("JK",names(peak_list),invert = TRUE)],
    function(x,all_peaks){
      ov <- findOverlaps(x,all_peaks)
      sqnms <- rep(NA,subjectLength(ov))
      sqnms[subjectHits(ov)] <- as.character(seqnames(x))
      dt <- data.table(seqnames = sqnms)
      dt[ subjectHits(ov), start := start(x) ]
      dt[ subjectHits(ov), end := end(x) ]
      return(dt)
    },all_peaks)


  out <- list(peaks = peaks, overlaps = overlaps , probs = probs)
  out <- c(out,peak_coordinates)

  return(out)
}

unified_lists <- unify_list(mosaics_gr)

save(unified_lists, file = "data/RData/unified_lists_wProbs.RData")

