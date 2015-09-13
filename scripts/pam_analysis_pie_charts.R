
rm(list = ls())

library(ggplot2)
library(reshape2)
library(data.table)

load( "data/RData/pam_analysis_K10.RData") ## results

dnase <- results[["dnase"]]

sets <- names(dnase)

states <- c("Active promoters","Strong enhancers","Weak enhancers",
  "Poised promoters","Heterochromatin","Others")
states <- data.table(id = 1:6,names = factor(states,levels = states))                       

remaining <- function(K,...)
{
  all <- list(...)
  if(class(all[[1]]) == "list"){
    all <- all[[1]]
  }
  all <- do.call(c,all)
  nop <- which(!1:K %in% all)
  if(length(nop) > 0){
    out <- (1:K)[nop]
  }else{
    out <- NULL
  }
  return(out)
}

pie_chart <- function(set,sizes,states)
{
  sizes <- data.table(clusters = names(sizes),sizes = sizes)
  setkey(sizes,clusters)
  nclust <- lapply(set,function(x,sizes)sizes[x,sum(sizes)],sizes)
  total <- sizes[,sum(sizes)]
  nclust <- do.call(c,nclust)
  out <- copy(states)
  out[,size := c(nclust, total - sum(nclust))]
  return(out)
}


## For EBNA3A:
## Active promoters - High H3K4me3 and H3K9ac are shown in clusters 1,5,7 and 9
## Strong enhancers - High H3K4me1 and H3K27ac are shown in clusters 2 and 8
## Weak enhancers - Intermediate H3K4me1 and Weak H3K27ac are shown in clusters 4 and 10
## Poised promoters - High H3K4me3 and low H3K27ac or high H3K27me3 which appear in cluster 3 and 10
## Heterochromatin - Cluster 6

ebna3a <- list()
ebna3a[[1]] <- c(1,5,7,9)
ebna3a[[2]] <- c(2,8)
ebna3a[[3]] <- c(4,10)
ebna3a[[4]] <- c(3)
ebna3a[[5]] <- 6

ebna3a <- pie_chart(ebna3a,results[["dnase"]]$EBNA3A$nr_clusters,states)

## For EBNA3B:
## Active promoters - High H3K4me3 and H3K9ac are shown in clusters 1 and 10
## Strong enhancers - High H3K4me1 and H3K27ac are shown in clusters 2,3 and 4
## Weak enhancers - Intermediate H3K4me1 and Weak H3K27ac are shown in clusters 5,7and 9
## Poised promoters - High H3K4me3 and low H3K27ac or high H3K27me3 which appear in cluster 6
## Heterochromatin - Cluster 8

ebna3b <- list()
ebna3b[[1]] <- c(1,10)
ebna3b[[2]] <- 2:4
ebna3b[[3]] <- c(5,7,9)
ebna3b[[4]] <- 6
ebna3b[[5]] <- 8

ebna3b <- pie_chart(ebna3b,results[["dnase"]]$EBNA3B$nr_clusters,states)

## For EBNA3C:
## Active promoters - High H3K4me3 and H3K9ac are shown in cluster 6
## Strong enhancers - High H3K4me1 and H3K27ac are shown in clusters 2, 4 and 9
## Weak enhancers - Intermediate H3K4me1 and Weak H3K27ac are shown in cluster 1
## Poised promoters - High H3K4me3 and low H3K27ac or high H3K27me3 which appear in cluster 3
## Heterochromatin - Cluster 7 and 8

ebna3c <- list()
ebna3c[[1]] <- 6
ebna3c[[2]] <- c(2,4,9)
ebna3c[[3]] <- 1
ebna3c[[4]] <- 3
ebna3c[[5]] <- 7:8

ebna3c <- pie_chart(ebna3c,results[["dnase"]]$EBNA3C$nr_clusters,states)

## For EBNA2:
## Active promoters - High H3K4me3 and H3K9ac are shown in clusters 2 and 5
## Strong enhancers - High H3K4me1 and H3K27ac are shown in clusters 6, 7 and 8
## Weak enhancers - Intermediate H3K4me1 and Weak H3K27ac are shown in clusters 1, 3 and 4
## Poised promoters - High H3K4me3 and low H3K27ac or high H3K27me3 which appear in cluster 10
## Heterochromatin - Cluster 9

ebna2 <- list()
ebna2[[1]] <- c(2,5)
ebna2[[2]] <- 6:8
ebna2[[3]] <- c(1,3,4)
ebna2[[4]] <- 10
ebna2[[5]] <- 9

ebna2 <- pie_chart(ebna2,results[["dnase"]]$EBNA2$nr_clusters,states)

## For RBPJ:
## Active promoters - High H3K4me3 and H3K9ac are shown in cluster 6
## Strong enhancers - High H3K4me1 and H3K27ac are shown in clusters 2,4,8 and 10
## Weak enhancers - Intermediate H3K4me1 and Weak H3K27ac are shown in clusters 1, 3 and 5
## Poised promoters - High H3K4me3 and low H3K27ac or high H3K27me3 which appear in cluster 7
## Heterochromatin - Cluster 9

rbpj <- list()
rbpj[[1]] <- 6
rbpj[[2]] <- c(2,4,8,10)
rbpj[[3]] <- c(1,3,5)
rbpj[[4]] <- 7
rbpj[[5]] <- 9

rbpj <- pie_chart(rbpj,results[["dnase"]]$RBPJ$nr_clusters,states)

base <- ggplot(ebna3b , aes(x = factor(1),y = size,fill = names),colour = "black")+
  geom_bar(stat = "identity",width = 1)+coord_polar(theta = "y")+
  theme_bw()+theme(axis.text = element_blank(),axis.ticks = element_blank(),
    panel.grid = element_blank())+scale_fill_brewer(palette = "Set1")+xlab("")+ylab("")
  

figs_dir <- "figures/for_paper"

pdf(file = file.path(figs_dir,"fig3b_ebna3b.pdf"))
print("EBNA3B")
print(ebna3b)
print(base + ggtitle("EBNA3B"))
dev.off()
    
pdf(file = file.path(figs_dir,"fig3b_ebna3a.pdf"))
print("EBNA3A")
print(ebna3a)
print(base + ggtitle("EBNA3A"))
dev.off()
    
pdf(file = file.path(figs_dir,"fig3b_ebna3c.pdf"))
print("EBNA3C")
print(ebna3c)
print(base + ggtitle("EBNA3C"))
dev.off()

pdf(file = file.path(figs_dir,"fig3b_ebna2.pdf"))
print("EBNA2")
print(ebna2)
print(base + ggtitle("EBNA2"))
dev.off()

pdf(file = file.path(figs_dir,"fig3b_rbpj.pdf"))
print("RBPJ")
print(rbpj)
print(base + ggtitle("RBPJ"))
dev.off()
    

