
rm(list = ls())

## library(Segvis)
library(parallel)

library(devtools)
load_all("~/Desktop/Docs/Code/Segvis")

load("data/ranges/all_EBV_GenomicRanges.RData")

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

fl <- 200   #### using the same 
base <- buildSegvis(name = "",file = "", maxBandwidth = 501,chr = "human",fragLen = fl)

segvis <- lapply(read_files[,(sets)],
                 function(x,base){
                   name(base) <- x
                   return(base)},base)

segvis <- mapply(function(seg,filename){
  file(seg) <- filename
  return(seg)},segvis,file.path(read_dir,read_files[,(files)]),SIMPLIFY = FALSE)

fix <- lapply(ranges,function(x){
  elementMetadata(x) <- NULL
  seqnames(x) <- as.character(seqnames(x))
  return(x)})

segvis <- mapply(function(seg,reg){
  regions(seg) <- reg
  return(seg)},segvis,fix[read_files[,(sets)]],SIMPLIFY = FALSE)
               
segvis <- lapply(segvis,loadReads,mc)

segvis <- lapply(segvis,matchReads,mc)
segvis <- lapply(segvis,getCoverage,mc)

summits <- lapply(segvis,findSummit,151,mc)

save(summits, file = "data/RData/summits.RData")

## to create rbpj set, we considered the peaks that are
## from jk234 first, and then added the ones in jk92 that didn't
## overlap jk234



    ## if(isPET(object)){
    ##   message("Setting PET flag")
    ##   ## when the reads are paired end tags, it gets the qname to match the
    ##   ## both ends of the fragment
    ##   pet_flag <- scanBamFlag(isPaired = TRUE)
    ##   param <- ScanBamParam(which = regions(object),flag = pet_flag,what = "qname")
    ## }else{
    ##   param <- ScanBamParam(which = regions(object))
    ## }    
    ## greads <- readGAlignments(file(object),
    ##   param = param,use.names = FALSE)    
    ## if(isPET(object)){
    ##   ## convert the qname into a numeric value for computation efficiency
    ##   qname <- as.numeric(as.factor(elementMetadata(greads)[["qname"]]))
    ##   greads <- .data.table.GRanges(as(greads, "GRanges"))
    ##   greads[,name:=qname] # add qname to greads
    ##   setorder(greads,name)
    ## }else{
    ##   greads <- .data.table.GRanges(as(greads, "GRanges"))
    ## }
    ## setkey(greads,seqnames,strand)    
    ## message("Bam file loaded")


## separate_reads <- function(x,chrom,str)x[seqnames == chrom & strand == str]

## bed2GRanges <- function(x){
##   gr <- GRanges(seqnames = x[,(seqnames)],
##                 ranges = IRanges( start = x[,(start)],
##                   end = x[,(end)]),strand = "*")
##   return(gr)
## }


## cross_corr <- function(chr,fwd,bwd,chrom.size,shift)
## {
##   message("...calculating strand cross correlation for ",chr)
##   fwd_cover <- coverage(trim(bed2GRanges(fwd[chr])),width = chrom.size[chr,(V2)])[[chr]]
##   bwd_cover <- coverage(trim(bed2GRanges(bwd[chr])),width = chrom.size[chr,(V2)])[[chr]]
##   out <- data.table(seqnames = chr,shift  = shift)
##   cc <- shiftApply(shift,fwd_cover,bwd_cover,cor,verbose =FALSE)
##   out[,strand.cross.corr := cc]
##   return(out)
## }

## dotprod <- function(strand.cross.corr,weights)return(sum(strand.cross.corr * weights))


## gr2data.table <- function(x)  data.table(seqnames = as.character(seqnames(x)),start = start(x),end = end(x))


## strand_cross_corr <- function(file,mc,chrom.size,shift = seq(50,300,by = 5))
## {
##   message("Processing ",file)
##   reads <- readGAlignments(file , param = NULL)
##   chr <- paste0("chr",c(1:19,"X","Y"))
##   message("Tidying data format")
##   reads <- as(reads,"GRanges") 
##   reads <- dropSeqlevels(reads,"chrM")
##   dt <- gr2data.table(reads)
##   dt[,strand := as.character(strand(reads))]
##   setkey(dt,strand)
##   fwd <- dt["+"]
##   bwd <- dt["-"]
##   setkey(fwd,seqnames)
##   setkey(bwd,seqnames)
##   fdepth <- fwd[,length(start),by = (seqnames)]
##   bdepth <- bwd[,length(start),by = (seqnames)]
##   if(length(chr) > nrow(fdepth)){
##     message(file, " contains less chromosomes than a complete genome, removing extra chromosome from chrom.size")
##     chr <- chr[chr %in% fdepth[,(seqnames)]]
##   }
##   rm(dt)
##   fwd[,end:=start]
##   bwd[,start:=end]
##   message("Calculating strand cross correlation")
##   dt <- mclapply(chr,cross_corr,fwd,bwd,chrom.size,shift = shift,mc.cores = mc,mc.preschedule=TRUE)
##   depth <- data.table(seqnames = fdepth[,(seqnames)], depth = fdepth[,(V1)] + bdepth[,(V1)])
##   dt <- do.call(rbind,dt)
##   rm(fwd,bwd,fdepth,bdepth)
##   depth[,weights:= depth / sum(depth)]
##   dt <- merge(dt,depth ,by = "seqnames")
##   du <- dt[,dotprod(strand.cross.corr,weights),by = (shift)]
##   rm(dt)
##   setnames(du,names(du),c("shift","strand.cross.corr"))
##   return(du)
## }

## cross_corr <- lapply(file.path(read_dir,read_files[,(files)]),
##   strand_cross_corr,mc,chrom.size,shift = seq(5,300,by = 5))
                     


