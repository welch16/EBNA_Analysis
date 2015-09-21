
rm(list = ls())

library(data.table)

load(file = "data/RData/unified_lists_wProbs.RData")
load(file = "data/RData/factor_overlaps.RData")

peaks <- unified_lists$peaks
overlaps <- unified_lists$overlaps
probs <- unified_lists$probs

out <- data.table(peaks,prob = probs$minProb,overlaps,factor_overlaps)

out_dir <- "inst/generated"
write.csv(file = file.path(out_dir,"supplement_matrix.csv"),out,row.names = FALSE,quote =FALSE)
