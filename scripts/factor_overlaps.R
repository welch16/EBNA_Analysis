rm(list = ls())

library(GenomicRanges)
library(data.table)
library(parallel)

load( file = "data/RData/unified_lists_wProbs.RData")

mc <- detectCores()

peaks <- unified_lists$peaks
overlaps <- unified_lists$overlaps
probs <- unified_lists$probs

sets <- c("EBNA2","EBNA3A","EBNA3B","EBNA3C","JK234","JK92")

tf_dir <- "inst/tf_peaks"
dhs_dir <- "inst/dnase_encode"

load_meta <- function(dir)
{
  files <- list.files(dir)
  meta <- files[grep("meta",files)]
  out <- data.table(read.table( file.path(dir,meta),sep = "\t",stringsAsFactors = FALSE,
                               header = TRUE))
  return(out)
}

tf_meta <- load_meta(tf_dir)
dhs_meta <- load_meta(dhs_dir)

load_files <- function(meta,dir,mc)
{  
  files <- paste0(meta[,(File.accession)],".bed.gz")
  bed_files <- mclapply(file.path(dir,files),read.table,
    stringsAsFactors = FALSE,header = FALSE,mc.cores = mc)
  gr <- lapply(bed_files,function(x){
    out <- GRanges(seqnames = x$V1,
                   ranges = IRanges(start = x$V2,end = x$V3),
                   strand = "*")
    return(out)})
  return(gr)
}

## only get optimal_idr peaks
tf_meta <- tf_meta[grep("optimal idr",Output.type)]
tf_peaks <- load_files(tf_meta,tf_dir,mc)
names(tf_peaks) <- tf_meta[,(File.accession)]

## only get peaks and not hotspots
dhs_meta <- dhs_meta[Output.type == "peaks"]
dhs_peaks <- load_files(dhs_meta,dhs_dir,mc)

dhs <- reduce(Reduce(c,dhs_peaks))

tfs <- tf_meta[,.(File.accession,Experiment.target)]
tfs[,Experiment.target:=gsub("-human","",Experiment.target)]
tfs <- split(tfs,tfs[,(Experiment.target)])

tfs <- mclapply(tfs,function(x,tf_peaks){
  peaks <- tf_peaks[x[,(File.accession)]]
  peaks <- reduce(Reduce(c,peaks))
  return(peaks)},tf_peaks,mc.cores = mc)
rm(tf_peaks)

gr <- GRanges(seqnames = peaks[,(seqnames)],
              ranges = IRanges(start = peaks[,(start)],
                end = peaks[,(end)]),
              strand = "*")
tfs[["Dnase"]] <- dhs

cols <- lapply(tfs,function(x,gr)
               ifelse(countOverlaps(gr,x) > 0,1,0),gr)

factor_overlaps <- as.data.table(cols)

save(factor_overlaps,file = "data/RData/factor_overlaps.RData")

