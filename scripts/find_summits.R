
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

sets <- c("EBNA2","EBNA3A","EBNA3B","EBNA3C","RBPJ","RBPJ")
segvis <- lapply(sets,
                 function(x,base){
                   name(base) <- x
                   return(base)},base)
names(segvis) <- read_files[,(sets)]

segvis <- mapply(function(seg,filename){
  file(seg) <- filename
  return(seg)},segvis,file.path(read_dir,read_files[,(files)]),SIMPLIFY = FALSE)


dt2gr <- function(x){
  GRanges(seqnames = x[,(seqnames)],
          ranges = IRanges(start = x[,(start)],end = x[,(end)]),
          strand = "*")}
            
fix <- lapply(sets,function(x,peaks){
  idx <- !is.na(peaks[[x]][,(seqnames)])
  out <- peaks[[x]][idx]
  out <- dt2gr(out)
  return(out)},unified_lists)

segvis <- mapply(function(seg,reg){
  regions(seg) <- reg
  return(seg)},segvis,fix,SIMPLIFY = FALSE)
               
segvis <- lapply(segvis,loadReads,mc)
segvis <- lapply(segvis,matchReads,mc)
segvis <- lapply(segvis,getCoverage,mc)


## using same bandwidth as in the histone profiless
ummits <- lapply(segvis,findSummit,151,mc)

summits <- lapply(read_files[,(sets)],function(x,summits,overlaps){
  if(x == "JK92" | x == "JK234"){
    out <- overlaps[["RBPJ"]]
  }else{
    out <- overlaps[[x]]
  }
  out[out == 1] <- summits[[x]]
  out[out == 0] <- NA
return(out)},summits,overlaps)

names(summits) <- paste0(read_files[,(sets)],".summit")
summits <- cbind(peaks[,1:3,with = FALSE],as.data.table(summits))

summits[,"RBPJ.summit" := ifelse(overlaps[["JK234"]] == 1,JK234.summit,JK92.summit)]

save(summits, file = "data/RData/summits.RData")

