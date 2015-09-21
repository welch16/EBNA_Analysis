
rm(list = ls())

library(data.table)
library(mosaics)
library(parallel)
library(GenomicRanges)

### set parameters
fragLen <- 200 
binSize <- 200 

read_dir <- "inst/bed_reads"
out_dir <- "inst/wig_files"
file_type <- "bed"

samples <- data.table(peaks = c("EBNA2",
                                "JK92",
                                "EBNA3C",
                                "JK234",
                                "EBNA3A",
                                "EBNA3B"),                        
                      files = c("s_1_92_seq_bowtie_uni",
                                "s_3_92_seq_bowtie_uni",
                                "s_6_234_seq_bowtie_uni",
                                "s_7_234_seq_bowtie_uni",
                                "s_1_250_seq_bowtie_uni",
                                "s_2_250_seq_bowtie_uni"))
setkey(samples,peaks)

samples[,files_fix := paste0(files,"_fix")]

chrom.sizes <- data.table(read.table("/p/keles/SOFTWARE/hg19.chrom.sizes"))

## re_create_bed <- function(set, samples,read_dir)
## {
##   dt <- data.table(read.table(file.path(read_dir,
##     paste0(samples[set,(files)],".bed")),stringsAsFactors = FALSE,header = FALSE,
##     comment.char = ""))
##   dt[,V4 := paste0("seq",1:nrow(dt))]
##   write.table(dt , file = file.path(read_dir,paste0(samples[set,(files_fix)],".bed")),
##     quote = FALSE, sep = "\t",col.names = FALSE,row.names = FALSE)
## }

##mclapply(samples[,(peaks)],re_create_bed,samples,read_dir,mc.cores = detectCores())

setkey(chrom.sizes,V1)

gen_wig_wrap <- function(set, samples, read_dir,out_dir,file_type, fl,bs,chrom.sizes)
{
  dt <- data.table(read.table(file = file.path(read_dir,
    paste0(samples[set,(files_fix)],".",file_type)),stringsAsFactors = FALSE,header =FALSE))

  gr <- GRanges(seqnames = dt[,(V1)],ranges = IRanges(start = dt[,(V2)],end = dt[,(V3)]),
                strand = dt[,(V6)])
  gr <- resize(gr,fl)
  bins <- lapply(as.character(chrom.sizes[,(V1)]),
        function(chr,bs,chrom.sizes){
          start <- seq(1,chrom.sizes[chr,(V2)] - bs,by = bs)         
          out <- IRanges(start , width = bs)
          return(GRanges(seqnames = chr, ranges = out ,strand = "*"))},
                 bs,chrom.sizes)
  bins <- do.call(c,bins)
  depth <- nrow(dt)
  mcols(bins)$wig <- round(1e6 *  countOverlaps(bins,gr) / depth,4)
  message("Start writing wig file for " , set )

  ff <- file.path(out_dir,paste0(set,".wig"))  
  for(chr in chrom.sizes[,(V1)]){
    message(chr)
    cat(paste0("variableStep chrom=",chr," span=",bs),file = ff ,append = TRUE , sep = "\n")
    bn <- subset(bins,seqnames == chr)
    lines <- paste(start(bn),bn$wig,sep = " ")
    z <- lapply(lines,cat , file = ff ,append = TRUE , sep = "\n")    
  }
}

lapply(samples[,(peaks)],gen_wig_wrap,samples,read_dir,out_dir,file_type,fragLen,binSize,chrom.sizes)

