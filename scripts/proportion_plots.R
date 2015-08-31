
rm(list = ls())

library(data.table)
library(ggplot2)

load("data/RData/unified_lists_wProbs.RData")
load("data/RData/factor_overlaps.RData")

all_data <- cbind(unified_lists[["overlaps"]],factor_overlaps)

sets <- c("EBNA2","EBNA3A","EBNA3B","EBNA3C","JK234","JK92","RBPJ")

get_proportions <- function(set,data,sets, DHS = FALSE)
{
  stopifnot(is.logical(DHS))
  idx <- which(data[[set]] == 1)
  my_data <- copy(data[idx])
  if(DHS){
    my_data <- my_data[Dnase == 1]
  }
  ## remove cols
  nms <- names(my_data)
  ## remove EBV peak-sets
  my_data <- my_data[ , ! nms %in% sets, with = FALSE]

  ## remove Dnase
  my_data[,Dnase := NULL]
  props <- colSums(my_data) / nrow(my_data)
  nms <- names(props)
  props <- data.table(set = set , TF = nms,prop = props)
  return(props)
}

props <- do.call(rbind,
  lapply(c("EBNA2","EBNA3A","EBNA3B","EBNA3C","RBPJ"),
         get_proportions,
         all_data,sets,FALSE))
props_dnase <- do.call(rbind,
  lapply(c("EBNA2","EBNA3A","EBNA3B","EBNA3C","RBPJ"),
         get_proportions,
         all_data,sets,TRUE))

proportions <- list(all = props, dnase = props_dnase)

save(proportions , file = "data/RData/TF_overlap_proportion.RData")


prop_plot <- function(set, data)
{
  idx <- data[["set"]] == set                
  data <- copy(data[idx])
  data <- data[order(-prop)]
  nms <- data[,(TF)]
  data[,TF := factor(TF,levels = nms)]
  p <- ggplot( data , aes(TF , prop, fill = TF))+ geom_bar(stat = "identity",width = .5,colour = "black",
                                       linesize = 1)+
    theme(axis.text.x = element_text(angle = 90, colour = "black",size = 6),legend.position = "none")+
    scale_fill_manual(values = rainbow(76))+
    xlab("ENCODE TFs")+ylab("Fraction of peaks")+ggtitle(set)+
    ylim(0,.85)
  return(p)
}


pdf(file = "figures/for_paper/fig6a.pdf")
print(prop_plot("EBNA3B",proportions[[1]]))
print(prop_plot("EBNA3B",proportions[[2]]))
dev.off()


pdf(file = "figures/for_paper/figS2_1.pdf")
print(prop_plot("EBNA3A",proportions[[1]]))
print(prop_plot("EBNA3A",proportions[[2]]))
dev.off()

pdf(file = "figures/for_paper/figS2_2.pdf")
print(prop_plot("EBNA3C",proportions[[1]]))
print(prop_plot("EBNA3C",proportions[[2]]))
dev.off()

pdf(file = "figures/for_paper/figS2_extra1.pdf")
print(prop_plot("EBNA2",proportions[[1]]))
print(prop_plot("EBNA2",proportions[[2]]))
dev.off()


pdf(file = "figures/for_paper/figS2_extra2.pdf")
print(prop_plot("RBPJ",proportions[[1]]))
print(prop_plot("RBPJ",proportions[[2]]))
dev.off()
