

rm( list = ls())

library(data.table)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(parallel)
library(pracma)

source("R/profile_functions.R")

profile_dir <- "data/RData/profile_matrices_old"

sets <- c("EBNA2","EBNA3A","EBNA3B","EBNA3C","RBPJ")

files <- list.files(profile_dir)
window_ext <- 2000
mc <- detectCores()

load("data/RData/unified_lists_wProbs.RData")
load("data/RData/factor_overlaps.RData")

profiles <- lapply(sets,
  function(x,dir){
    files <- grepV(x,files)
    out <- mclapply(file.path(dir,files),
      load_clean_profile,window_ext,mc.cores = mc)
    names(out) <- files
    return(out)},profile_dir)
names(profiles) <- sets

peak_list <- lapply(sets,add_match,unified_lists,factor_overlaps)
names(peak_list) <- sets


M <- 200
lb <- -M
ub <- M

stats <- lapply(profiles,function(x,lb,ub){
  out <- mclapply(x,function(y,lb,ub){
    y[between(coord,lb,ub),max(tagCounts),by  = match]
  },lb,ub,mc.cores = mc)
  return(out)},lb,ub)

stats <- lapply(profiles,function(x,lb,ub){
  out <- mclapply(x,function(y,lb,ub){
    y[between(coord,lb,ub),trapz(coord,tagCounts),by  = match]
  },lb,ub,mc.cores = mc)
  return(out)},lb,ub)

stats <- lapply(sets,add_names,stats)
names(stats) <- sets

stats <- do.call(rbind,stats)

mat <- peak_list[["EBNA3A"]][Dnase == 1][,(match)]

ggplot(stats[match %in% mat & rep <= 2],aes(histone  , V1,color = as.factor(rep)))+geom_boxplot()+
  scale_y_log10()+facet_grid( set ~ .)+scale_color_brewer(palette = "Set1")+
  theme(axis.text.x = element_text(angle = 90))
dev.off()




setkey(stats,set)
test <- copy(stats["EBNA3A"][match %in% mat & rep == 1])
test <- dcast.data.table(data =test, match ~ histone,value.var = "V1") 

heatmap(as.matrix(test[,-1,with = FALSE]),col = brewer.pal(10,"RdYlGn"),labRow = "")
dev.off()





## test <- list()
## suppressWarnings(
## test[[set]] <- mclapply(profiles[["EBNA3B"]],
##   function(x){    
##     out <- x[,mean(tagCounts,trim = .05,na.rm = TRUE),by = coord]
##     return(out)},mc.cores = mc)
##                  )
## test <- add_names("EBNA3B",test,files)
## test <- do.call(rbind,test)


## pdf(ile = "EBNA3B_old.pdf")
## ggplot(test[grep("H3",histone)],
##        aes(coord,V1,colour = histone,linetype = as.factor(rep)))+
##   geom_line(linesize = 1)+scale_color_brewer(palette = "Set1")+
##   theme(panel.background = element_rect(fill = "black"),
##         panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ylim(0,2)
## dev.off()


