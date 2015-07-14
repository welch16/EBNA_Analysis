

rm(list = ls())

library(GenomicRanges)
library(data.table)
library(parallel)
library(reshape2)
library(ggplot2)

load("data/ranges/all_EBV_GenomicRanges.RData") ## ranges

sets <- c("EBNA2","EBNA3A","EBNA3B","EBNA3C","RBPJ")


peakSets <- lapply(ranges,reduce)
peakSets <- peakSets[names(peakSets) %in% sets]


bigSet <- peakSets[[1]]
for(k in 2:length(peakSets)){
  bigSet <- union(bigSet,peakSets[[k]])
}

bigSet <- reduce(bigSet)

mc <- detectCores()

incidence <- mclapply(peakSets,function(x,big){
  ifelse(countOverlaps(big,x) >0,1,0)},bigSet,mc.cores = mc)

elementMetadata(bigSet) <- incidence

extra <- names(elementMetadata(ranges[[1]]))
## remove annot columns (no way to re-estimate those for the new regions)
extra <- extra[-(1:4)]


extra <- extra[!extra %in% sets]

## remove allTF
extra <- extra[-length(extra)]

overlaps <- mclapply(peakSets,findOverlaps,bigSet,mc.cores = mc)

get_extra_col <- function(name,overlaps,peakSets,bigSet)
{
  ## extract extra cols
  extra_cols <- lapply(peakSets,function(x,nm)elementMetadata(x)[[nm]],name)
  extra_cols <- lapply(extra_cols,as.numeric)
  extra_cols <- mapply(function(extra,ov,n){
    out <- rep(0,n)
    out[subjectHits(ov)] <- extra
    return(out)
  },extra_cols,MoreArgs= list(length(bigSet)),overlaps,SIMPLIFY=FALSE)
  out <- rep(0,length(bigSet))
  for(k in 1:length(extra_cols)){
    out <- out + extra_cols[[k]]
  }
  out <- ifelse(out > 0,1,0)
  return(out)
}

extra_cols <- mclapply(extra,get_extra_col,overlaps,
       ranges[names(ranges) %in% sets],bigSet,mc.cores =mc)
names(extra_cols) <- extra

extra_cols <- lapply(extra_cols,as.matrix)
extra_mat<- do.call(cbind,extra_cols)

dt <- data.table(as.data.frame(bigSet))

colnames(extra_mat) <- extra
dt <- cbind(dt,extra_mat)

out_dir <- "inst/generated"
write.csv(file = file.path(out_dir,"Big_overlaps_matrix.csv"),dt,row.names = FALSE,quote =FALSE)
                   
