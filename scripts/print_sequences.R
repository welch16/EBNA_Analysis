
rm(list = ls())

library(data.table)

load("data/RData/sequences_around_summit.RData")
load("data/RData/unified_lists_wProbs.RData")
load("data/RData/factor_overlaps.RData")

peaks <- unified_lists$peaks
overlaps <- unified_lists$overlaps
probs <- unified_lists$probs

sets <- c("EBNA2","EBNA3A","EBNA3B","EBNA3C","RBPJ")

prepare_data <- function(set, peaks,overlaps,probs,factor_overlaps,sequences)
{
  idx <- which(overlaps[[set]] == 1)
  out <- peaks[idx]
  out[[paste0(set,".prob")]] <- probs[[set]][idx]
  out[["minProb"]] <- probs[["minProb"]][idx]
  out[["Dnase"]] <- factor_overlaps[["Dnase"]][idx]
  seq <- sequences[[paste0(set,".sequence")]]
  seq <- seq[!is.na(seq)]
  out[["sequence"]] <- seq
  return(out)
}

data_sets <- lapply(sets,prepare_data,peaks,overlaps,probs,factor_overlaps,sequences)

write <- function(set,data)
{
  out_dir <- "data/txt_files/"
  write.table(data , file = file.path(out_dir,paste0(set,"_sequences.txt")),
              sep = "\t",col.names = TRUE,quote = FALSE,row.names = FALSE)
}

mapply(write , sets, data_sets,SIMPLIFY = FALSE)
