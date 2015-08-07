
rm(list = ls())

library(GenomicAlignments)
library(GenomicRanges)
library(data.table)
library(dplyr)
library(ggplot2)
library(reshape2)
library(parallel)
library(rtracklayer)
library(RColorBrewer)

mc <- detectCores()

load("data/ranges/all_EBV_GenomicRanges.RData")

chromHMM <- read.table("inst/chromHMM/wgEncodeBroadHmmGm12878HMM.bed",header = FALSE,
                       stringsAsFactors = FALSE)
chromHMM <- data.table(chromHMM[,-(5:8)])

extra <- chromHMM[,(V9)]
chromHMM <- chromHMM[,1:4,with  = FALSE]
setnames(chromHMM,names(chromHMM),c("seqnames","start","end","label"))

chromHMM_gr <- GRanges(seqnames = chromHMM[,(seqnames)],
                       ranges = IRanges(
                         start = chromHMM[,(start)],
                         end = chromHMM[,(end)]),
                       strand = "*")
chromHMM_gr$label <- chromHMM[,(label)]
chain <- import.chain("inst/chromHMM/hg18ToHg19.over.chain")

hg19_annot <- liftOver(chromHMM_gr,chain)
hg19_annot <- unlist(hg19_annot)

ranges <- lapply(ranges,subset , Dnase == 1)

overlaps <- mclapply(ranges,findOverlaps,hg19_annot,mc.cores = mc)


hg19_dt <- data.table(seqnames = as.character(seqnames(hg19_annot)),
                      start = start(hg19_annot),
                      end = end(hg19_annot),
                      label = hg19_annot$label)

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

chromHMM_label <- function(ebv_set,ov,hg19_dt,rule)
{
  ebv_label <- ebv_set[queryHits(ov)]
  ebv_label <- data.table(seqnames = as.character(seqnames(ebv_label)),
                          start = start(ebv_label),
                          end = end(ebv_label))
  ebv_label[,id := paste0(seqnames,":",start,"-",end)]
  to_add <- copy(hg19_dt)
  setnames(to_add,names(to_add),paste0("chromHMM_",names(to_add)))
  ebv_label <- cbind(ebv_label,to_add[subjectHits(ov)])
  ebv_label[,chromHMM_width := chromHMM_end - chromHMM_start + 1]
  ext <- ebv_label[ , untie(chromHMM_label,chromHMM_width , rule),by = id]
  ebv <- data.table(seqnames = as.character(seqnames(ebv_set)),                    
                    start = start(ebv_set),end = end(ebv_set))
  ebv[,id := paste0(seqnames,":",start,"-",end)]
  ebv <- merge(ebv,ext,by = "id",all.x = TRUE)
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
                   MoreArgs = list(hg19_dt = hg19_dt, rule = "maxwidth"),
                   SIMPLIFY = FALSE,mc.cores = mc)

dt <- mapply(function(x,y){
  data.table(set = x,label = y)},names(labels_mw),labels_mw,SIMPLIFY = FALSE)


labs <- c("0_none",
          "1_Active_Promoter",
          "2_Weak_Promoter",
          "3_Poised_Promoter",
          "4_Strong_Enhancer",
          "5_Strong_Enhancer",
          "6_Weak_Enhancer",
          "7_Weak_Enhancer",
          "8_Insulator",
          "9_Txn_Transition",
          "10_Txn_Elongation",
          "11_Weak_Txn",
          "12_Repressed",
          "13_Heterochrom/lo",
          "14_Repetitive/CNV",
          "15_Repetitive/CNV")

dt <- do.call(rbind, dt)
dt[is.na(label),label := "0_none"]

dt[,label := factor(label,levels = labs)]

new_labs <- c("0_none",
          "1_Active_Promoter",
          "2_Weak_Promoter",
          "3_Poised_Promoter",
          "4_Strong_Enhancer",
          "4_Strong_Enhancer",
          "6_Weak_Enhancer",
          "6_Weak_Enhancer",
          "8_Insulator",
          "9_Txn_Transition",
          "10_Txn_Elongation",
          "11_Weak_Txn",
          "12_Repressed",
          "13_Heterochrom/lo",
          "14_Repetitive/CNV",
          "14_Repetitive/CNV")


dt[,label := plyr::mapvalues(label ,
      from = labs, to = new_labs)]

dt[,label := plyr::revalue(label,
  c("6_Weak_Enhancer"="5_Weak_Enhancer",
    "8_Insulator"="6_Insulator",
    "9_Txn_Transition"="7_Txn_Transition",
    "10_Txn_Elongation"="8_Txn_Elongation",
    "11_Weak_Txn"="9_Weak_Txn",
    "12_Repressed"="10_Repressed",
    "13_Heterochrom/lo"="11_Heterochrom/lo",
    "14_Repetitive/CNV"="12_Repetitive/CNV"))]      

tab <- table(dt)
du <- data.table(tab)

labs <- unique(du[,(label)])

du[,label := factor(label, levels = rev(labs))]
du[,perc := 100 * N / sum(N) , by  = set]


getPalette1 <- colorRampPalette(brewer.pal(8,"Set1"))
getPalette2 <- colorRampPalette(brewer.pal(8,"Set2"))
vals <- c( getPalette1(6), getPalette2(7))

p <- ggplot(du[grep("JK",set,invert = TRUE)] , aes(label,perc,fill = label ))+geom_bar(width = .8,stat = "identity")+
  scale_fill_manual(values = rev(vals))+scale_y_continuous(expand = c(0,0),limits = c(0,70))+
  coord_flip()+facet_grid( set ~ .)+theme(legend.position = "none")+ylab("Percentage")+
  xlab("ChromHMM annotation")


sets <- names(ranges)
sets <- sets[grep("JK",sets,invert = TRUE)]

du[,label := factor(label, levels = labs)]

pies <- lapply( sets , function(z){
  ggplot(du[set == z] , aes(x = factor(1),y = perc, fill = label),colour = "black")+
    geom_bar(stat = "identity",width = 1) + coord_polar(theta = "y")+
    scale_fill_manual(name="Category",values = vals)+
    theme_bw()+theme(axis.text = element_blank(),axis.ticks = element_blank(),
                     panel.grid = element_blank())+
    xlab("")+ylab("")})
names(pies) <- sets

pies <- mapply(function(x,y){
  x + ggtitle(y)},pies,sets,SIMPLIFY = FALSE)


chromHMM <- list(table = tab, data = du , plot = p,pies = pies)

save(file = "data/RData/chromHMM_proportions_dnase.RData",chromHMM)
