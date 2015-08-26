
rm(list = ls())

library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(data.table)
library(parallel)

load( file = "data/RData/summits.RData")
load( file = "data/RData/unified_lists_wProbs.RData")

seq_ext <- 250
mc <- detectCores()

peaks <- unified_lists$peaks
overlaps <- unified_lists$overlaps
probs <- unified_lists$probs

sets <- c("EBNA2","EBNA3A","EBNA3B","EBNA3C","JK234","JK92")

chrom.sizes <-read.table("/p/keles/SOFTWARE/hg19.chrom.sizes",
                         stringsAsFactors = FALSE)
chrom.sizes <- data.table(chrom.sizes)

setnames(chrom.sizes,names(chrom.sizes),c("seqnames","seqlength"))
setkey(chrom.sizes,seqnames)
chrom.sizes <- chrom.sizes[!"chrY"]


extract_sequences <- function(set,summits,chrom.sizes,seq_ext)
{
  col <- paste0(set,".summit")
  summ <- summits[[col]]
  idx <- which(!is.na(summ))  
  regions <- GRanges( seqnames = summits[idx,(seqnames)],
                     ranges = IRanges(start = summ[idx] - seq_ext,
                       end = summ[idx] + seq_ext),
                     strand = "*")
  regions <- split(regions, seqnames(regions))
  sequences <- mclapply( chrom.sizes[,(seqnames)],
    function(chr,regions,chrom.sizes){
      reg <- regions[[chr]]
      len <- chrom.sizes[chr,(seqlength)]
      idx <- which(end(reg) < len)
      seq <- getSeq(  Hsapiens,reg[idx],as.character = TRUE)
      out <- rep("",length(reg))
      out[idx] <- seq
      out[-idx] <- do.call(paste0,as.list(rep("N",seq_ext * 2 + 1)))
      return(out)},regions,chrom.sizes,mc.cores = mc)
  names(sequences) <- chrom.sizes[,(seqnames)]
  out <- copy(summits[,c("seqnames","start","end",col),with = FALSE])
  nms <- gsub("summit","sequence",names(out))
  setnames(out,names(out),nms)
  setkey(out,seqnames)
  col <- paste0(set,".sequence")
  seq <- lapply(chrom.sizes[,(seqnames)],
                function(ch,sequences,out){
                  output <- out[ch][[col]]
                  idx <- which(!is.na(output))
                  output[idx] <- sequences[[ch]]
                  return(output)},sequences,out)
  names(seq) <- chrom.sizes[,(seqnames)]
  out[[col]] <- ""
  for(ch in chrom.sizes[,(seqnames)]){
    out[ch][[col]] <- seq[[ch]]
  }  
  return(out)
}

sequences <- lapply(sets,extract_sequences,summits,chrom.sizes,seq_ext)
sequences <- do.call(cbind,lapply(sequences,function(x)x[,4,with = FALSE]))

sequences[,RBPJ.sequence := ifelse(is.na(JK234.sequence), JK92.sequence,JK234.sequence)]

save(sequences,file = "data/RData/sequences_around_summit.RData")
