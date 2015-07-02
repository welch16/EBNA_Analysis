
rm(list = ls())

library(GenomicRanges)
library(data.table)
library(broom)

chrom.info <- read.table(file="/p/keles/SOFTWARE/hg19.chrom.sizes", header=FALSE)
names(chrom.info) <- c('chrom','length')
chrom.info$is_circular <- rep(FALSE, dim(chrom.info)[1])
chrom.info <- chrom.info[-21,]

load("data/ranges/all_EBV_GenomicRanges.RData")

re_annotate <- function(col)
{
  out <- rep(0,3)
  out[1] <- col[["5p1"]]
  out[2] <- col[["gene"]]
  out[3] <- sum(col) - sum(out)
  names(out) <- c("Promoters","Gene body","Intergenic")
  return(out)
}


association_test <- function(ranges,peak_set1,peak_set2){
  col1 <- table(ranges[[peak_set1]]$region)
  col2 <- table(ranges[[peak_set2]]$region)  
  cont_tab<- rbind(re_annotate(col1),re_annotate(col2))
  rownames(cont_tab) <- c(peak_set1,peak_set2)
  out <- c(Set1=peak_set1,Set2=peak_set2)
  out <- data.table(cbind(t(out),tidy(chisq.test(cont_tab))))
  return(out)  
}

k <- 1
assoc <- list()
peak_sets <- names(ranges)
for(j in 1:(length(ranges)-1)){
  for(i in (j+1):length(ranges)){
    assoc[[k]] <- association_test(ranges,peak_sets[j],peak_sets[i])
    k <- k+1
  }
}

assoc <- do.call(rbind,assoc)
assoc[,minus_log10 := -log10(p.value)]

save(assoc,  file = "data/RData/association_tests.RData" )
