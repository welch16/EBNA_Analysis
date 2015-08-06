
rm(list = ls())

library(GenomicAlignments)
library(GenomicRanges)
library(data.table)
library(dplyr)
library(ggplot2)
library(reshape2)
library(parallel)

mc <- detectCores()


load("data/ranges/all_EBV_GenomicRanges.RData")

chromHMM <- read.table("inst/chromHMM/GM12878_chromHMM_annot.bed",header = FALSE)
chromHMM <- data.table(chromHMM[,-(5:8)])

extra <- chromHMM[,(V9)]
chromHMM <- chromHMM[,1:4,with  = FALSE]
setnames(chromHMM,names(chromHMM),c("seqnames","start","end","label"))

chromHMM_gr <- GRanges(seqnames = chromHMM[,(seqnames)],
                       ranges = IRanges(
                         start = chromHMM[,(start)],
                         end = chromHMM[,(end)]),
                       strand = "*")

ranges <- lapply(ranges,subset, Dnase == 1)

overlaps <- mclapply(ranges,findOverlaps,chromHMM_gr,mc.cores = mc)

untie <- function(chromHMM_label,chromHMM_width, rule)
{
  if(length(chromHMM_label) == 1){
    out <- chromHMM_label
  }else{
    stopifnot(rule %in% c("maxwidth","first","last"))
    if(rule == "maxwidth"){
      out <- chromHMM_label[which.max(chromHMM_width)]
    }else if(rule == "first"){     
      out <- chromHMM_label[1]     
    }else if(rule == "last"){
      out <- chromHMM_label[length(chromHMM_label)]
    }
  }
  return(out)
}

chromHMM_label <- function(ebv_set,ov,chromHMM,rule)
{
  ebv_label <- ebv_set[queryHits(ov)]
  ebv_label <- data.table(seqnames = as.character(seqnames(ebv_label)),
                          start = start(ebv_label),
                          end = end(ebv_label))
  ebv_label[,id := paste0(seqnames,":",start,"-",end)]
  to_add <- copy(chromHMM)
  setnames(to_add,names(to_add),paste0("chromHMM_",names(to_add)))
  ebv_label <- cbind(ebv_label,to_add[subjectHits(ov)])
  ebv_label[,chromHMM_width := chromHMM_end - chromHMM_start + 1]
  ext <- ebv_label[ , untie(chromHMM_label,chromHMM_width , rule),by = id]
  ebv <- data.table(seqnames = as.character(seqnames(ebv_set)),                    
                    start = start(ebv_set),end = end(ebv_set))
  ebv[,id := paste0(seqnames,":",start,"-",end)]
  ebv <- merge(ebv,ext,by = "id")
  ebv_set_new <- GRanges(seqnames = ebv[,(seqnames)],
                         ranges = IRanges(
                           start = ebv[,(start)],
                           end = ebv[,(end)]),
                         strand = "*")
  ebv_set_new$chromHMM <- ebv[,(V1)]
  ov <- findOverlaps(ebv_set,ebv_set_new)
  out <- ebv_set_new$chromHMM[subjectHits(ov)]
  return(out)                        
}


labels_mw <- mcmapply(chromHMM_label,
                   ranges,overlaps,
                   MoreArgs = list(chromHMM = chromHMM, rule = "maxwidth"),
                   SIMPLIFY = FALSE,mc.cores = mc)

dt <- mapply(function(x,y){
  data.table(set = x,label = y)},names(labels_mw),labels_mw,SIMPLIFY = FALSE)

dt <- do.call(rbind, dt)
dt[,label := factor(label)]

labs <- c("TSS","PF","E","WE","CTCF","T","R")
dt[,label := plyr::mapvalues(label ,
      from = c("1","2","3","4","5","6","7"),
      to = labs)]

tab <- table(dt)

du <- data.table(tab)
du[,label := factor(label, levels = rev(labs))]

du[,perc := 100 * N / sum(N) , by  = set]

p <- ggplot(du[grep("JK",set,invert = TRUE)] , aes(label,perc,fill = label ))+geom_bar(width = .8,stat = "identity")+
  scale_fill_brewer(palette = "Set1")+scale_y_continuous(expand = c(0,0),limits = c(0,50))+
  coord_flip()+facet_grid( set ~ .)+theme(legend.position = "none")+ylab("Percentage")+
  xlab("ChromHMM annotation")

chromHMM <- list(table = tab, data = du , plot = p)

save(file = "data/RData/chromHMM_proportions_dnase.RData",chromHMM)
