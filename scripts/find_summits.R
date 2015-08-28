
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

peaks <- unified_lists$peaks
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

gr2dt <- function(x){
  dt <- data.table(seqnames = as.character(seqnames(x)),
             start = start(x),
             end = end(x))
  return(dt)
}

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
summits <- lapply(segvis,findSummit,151,mc)

summits <- lapply(read_files[,(sets)],function(x,summits,peaks,overlaps,segvis){

  if(x == "JK92" | x == "JK234"){
    x1 <- "RBPJ"
  }else{
    x1 <- x
  }

  reg <- gr2dt(regions(segvis[[x]]))
  reg[,summit := summits[[x]]]

  idx <-which( overlaps[[x1]] == 1)

  ov <- findOverlaps(dt2gr(peaks[idx]),dt2gr(reg))

  reordered_summit <- reg[subjectHits(ov) , (summit)]

  out <- overlaps[[x1]]
  out[out == 1] <- reordered_summit
  out[out == 0] <- NA
return(out)},summits,peaks,overlaps,segvis)

names(summits) <- paste0(read_files[,(sets)],".summit")
summits <- cbind(peaks[,1:3,with = FALSE],as.data.table(summits))

summits[,"RBPJ.summit" := ifelse(overlaps[["JK234"]] == 1,JK234.summit,JK92.summit)]

save(summits, file = "data/RData/summits.RData")

