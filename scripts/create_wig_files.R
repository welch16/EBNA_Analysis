
rm(list = ls())

library(data.table)
library(mosaics)


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

gen_wig_wrap <- function(set, samples, read_dir,out_dir,file_type, fl,bs)
{  
  generateWig(infile = file.path(read_dir,paste0(samples[set,(files)],".",file_type)),
              fileFormat = file_type,
              outfileLoc = out_dir,PET = FALSE,fragLen = fl,span = bs)             
}

lapply(samples$peaks,gen_wig_wrap,samples,read_dir,out_dir,file_type,fragLen,binSize)

