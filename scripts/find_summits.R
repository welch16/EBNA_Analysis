
rm(list = ls())

library(Segvis)
library(parallel)

## library(devtools)
## load_all("~/Desktop/Docs/Code/Segvis")

load( file = "data/RData/unified_lists_wProbs.RData")


read_files <- data.table(
     sets = c(
       "EBNA2",
       "EBNA3A",
       "EBNA3B",
       "EBNA3C",
       "JK92",
       "JK234"),
     files = c(
       "s_1_92_seq_bowtie_uni.bam",
       "s_1_250_seq_bowtie_uni.bam",
       "s_2_250_seq_bowtie_uni.bam",
       "s_6_234_seq_bowtie_uni.bam",
       "s_3_92_seq_bowtie_uni.bam",
       "s_7_234_seq_bowtie_uni.bam"))

read_dir <- "inst/reads"
setkey(read_files,sets)

chrom.size <- data.table(read.table("/p/keles/EnhancerPred/volumeC/mm10_chrom.sizes"))
setkey(chrom.size,V1)

mc <- detectCores()

peaks <- cbind(unified_lists$peaks,unified_lists$overlaps)
overlaps <- unified_lists$overlaps

fl <- 200   #### using the same fragment length used to call peaks
base <- buildSegvis(name = "",file = "", maxBandwidth = 501,chr = "human",fragLen = fl)

segvis <- lapply(read_files[,(sets)],
                 function(x,base){
                   name(base) <- x
                   return(base)},base)

segvis <- mapply(function(seg,filename){
  file(seg) <- filename
  return(seg)},segvis,file.path(read_dir,read_files[,(files)]),SIMPLIFY = FALSE)
names(segvis) <- read_files[,(sets)]

fix <- lapply(segvis,function(x,peaks){
  nm <- name(x)
  idx <- peaks[[nm]] == 1
  out <- peaks[idx]
  return(out)},peaks)

segvis <- mapply(function(seg,reg){
  regions(seg) <- GRanges(seqnames = reg[,(seqnames)],
                          ranges = IRanges(start = reg[,(start)],end = reg[,(end)]),
                          strand = "*")
  return(seg)},segvis,fix,SIMPLIFY = FALSE)
               
segvis <- lapply(segvis,loadReads,mc)
segvis <- lapply(segvis,matchReads,mc)
segvis <- lapply(segvis,getCoverage,mc)


## using same bandwidth as in the histone profiles
summits <- lapply(segvis,findSummit,151,mc)

summits <- lapply(read_files[,(sets)],function(x,summits,overlaps){
  out <- overlaps[[x]]
  out[out == 1] <- summits[[x]]
  out[out == 0] <- NA
  return(out)},summits,overlaps)
names(summits) <- paste0(read_files[,(sets)],".summit")

summits <- cbind(peaks[,1:3,with = FALSE],as.data.table(summits))

summits[,"RBPJ.summit" := ifelse(!is.na(JK234.summit),JK234.summit,JK92.summit)]

save(summits, file = "data/RData/summits.RData")

