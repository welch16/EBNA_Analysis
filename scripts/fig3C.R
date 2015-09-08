

rm( list = ls())

library(data.table)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(parallel)
library(pracma)
library(cluster)

source("R/profile_functions.R")

profile_dir <- "data/RData/profile_matrices_old"

sets <- c("EBNA2","EBNA3A","EBNA3B","EBNA3C","RBPJ")
sets <- c("EBNA3A","EBNA3B","EBNA3C")


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

profiles <- lapply(sets,add_names,profiles)

profiles <- lapply(profiles,function(x)x[rep == 1])

profiles <- do.call(rbind,profiles)

set_dnase <- lapply(sets,
  function(set,peak_list)peak_list[[set]][Dnase == 1,(match)],
    peak_list)
names(set_dnase) <- sets

profiles_dnase <- lapply(sets,
  function(ss,set_dnase,stats){
    out <- stats[set == ss]
    matchs <- set_dnase[[ss]]
    setkey(out,match)
    return(out[matchs])},set_dnase,profiles)
profiles_dnase <- do.call(rbind,profiles_dnase)

setkey(profiles,histone)
setkey(profiles_dnase,histone)

histones <- profiles[,unique(histone)]
histones <- grepV("H3",histones)
histones <- grepV("79",histones,invert = TRUE)

profiles <- profiles[histones]
profiles_dnase <- profiles_dnase[histones]

profiles[, rep := NULL]
profiles_dnase[,rep := NULL]

trim <- .05
suppressWarnings(
mean_profiles <- profiles[,mean(tagCounts,trim,na.rm = TRUE),by = .(coord,set,histone)]
                )

suppressWarnings(
mean_profiles_dnase <- profiles_dnase[,mean(tagCounts,trim,na.rm = TRUE),by = .(coord,set,histone)]
                 )

figs_dir <- "figures/for_paper"

setkey(mean_profiles,set)
setkey(mean_profiles_dnase,set)

pdf(file = file.path(figs_dir, "fig3C_EBNA3B.pdf")  ,width = 9,height = 8)
ggplot(mean_profiles["EBNA3B"],aes(coord,V1,colour = histone))+geom_line(size = 1.2)+
  scale_color_brewer(palette = "Dark2")+
  theme(legend.position = "bottom")+geom_vline(xintercept = 0 , linetype =2 ,size = 1.2)+
  guides(colour = guide_legend(nrow = 2))+ylim(0,2)+ylab("Normalized signal")+
  xlab("Distance to peaks summit")
ggplot(mean_profiles_dnase["EBNA3B"],aes(coord,V1,colour = histone))+geom_line(size = 1.2)+
  scale_color_brewer(palette = "Dark2")+
  theme(legend.position = "bottom")+geom_vline(xintercept = 0 , linetype =2 ,size = 1.2)+
  guides(colour = guide_legend(nrow = 2))+ylim(0,2)+ylab("Normalized signal")+
  xlab("Distance to peaks summit")
dev.off()


pdf(file = file.path(figs_dir, "fig3C_EBNA3A.pdf")  ,width = 9,height = 8)
ggplot(mean_profiles["EBNA3A"],aes(coord,V1,colour = histone))+geom_line(size = 1.2)+
  scale_color_brewer(palette = "Dark2")+
  theme(legend.position = "bottom")+geom_vline(xintercept = 0 , linetype =2 ,size = 1.2)+
  guides(colour = guide_legend(nrow = 2))+ylim(0,2)+ylab("Normalized signal")+
  xlab("Distance to peaks summit")
ggplot(mean_profiles_dnase["EBNA3A"],aes(coord,V1,colour = histone))+geom_line(size = 1.2)+
  scale_color_brewer(palette = "Dark2")+
  theme(legend.position = "bottom")+geom_vline(xintercept = 0 , linetype =2 ,size = 1.2)+
  guides(colour = guide_legend(nrow = 2))+ylim(0,2)+ylab("Normalized signal")+
  xlab("Distance to peaks summit")
dev.off()


pdf(file = file.path(figs_dir, "fig3C_EBNA3C.pdf")  ,width = 9,height = 8)
ggplot(mean_profiles["EBNA3C"],aes(coord,V1,colour = histone))+geom_line(size = 1.2)+
  scale_color_brewer(palette = "Dark2")+
  theme(legend.position = "bottom")+geom_vline(xintercept = 0 , linetype =2 ,size = 1.2)+
  guides(colour = guide_legend(nrow = 2))+ylim(0,2)+ylab("Normalized signal")+
  xlab("Distance to peaks summit")
ggplot(mean_profiles_dnase["EBNA3C"],aes(coord,V1,colour = histone))+geom_line(size = 1.2)+
  scale_color_brewer(palette = "Dark2")+
  theme(legend.position = "bottom")+geom_vline(xintercept = 0 , linetype =2 ,size = 1.2)+
  guides(colour = guide_legend(nrow = 2))+ylim(0,2)+ylab("Normalized signal")+
  xlab("Distance to peaks summit")
dev.off()



## pdf(file = file.path(figs_dir, "fig3C_EBNA2.pdf")  ,width = 9,height = 8)
## ggplot(mean_profiles["EBNA2"],aes(coord,V1,colour = histone))+geom_line(size = 1.2)+
##   scale_color_brewer(palette = "Dark2")+
##   theme(legend.position = "bottom")+geom_vline(xintercept = 0 , linetype =2 ,size = 1.2)+
##   guides(colour = guide_legend(nrow = 2))+ylim(0,2)+ylab("Normalized signal")+
##   xlab("Distance to peaks summit")
## ggplot(mean_profiles_dnase["EBNA2"],aes(coord,V1,colour = histone))+geom_line(size = 1.2)+
##   scale_color_brewer(palette = "Dark2")+
##   theme(legend.position = "bottom")+geom_vline(xintercept = 0 , linetype =2 ,size = 1.2)+
##   guides(colour = guide_legend(nrow = 2))+ylim(0,2)+ylab("Normalized signal")+
##   xlab("Distance to peaks summit")
## dev.off()

## pdf(file = file.path(figs_dir, "fig3C_RBPJ.pdf")  ,width = 9,height = 8)
## ggplot(mean_profiles["RBPJ"],aes(coord,V1,colour = histone))+geom_line(size = 1.2)+
##   scale_color_brewer(palette = "Dark2")+
##   theme(legend.position = "bottom")+geom_vline(xintercept = 0 , linetype =2 ,size = 1.2)+
##   guides(colour = guide_legend(nrow = 2))+ylim(0,2)+ylab("Normalized signal")+
##   xlab("Distance to peaks summit")
## ggplot(mean_profiles_dnase["RBPJ"],aes(coord,V1,colour = histone))+geom_line(size = 1.2)+
##   scale_color_brewer(palette = "Dark2")+
##   theme(legend.position = "bottom")+geom_vline(xintercept = 0 , linetype =2 ,size = 1.2)+
##   guides(colour = guide_legend(nrow = 2))+ylim(0,2)+ylab("Normalized signal")+
##   xlab("Distance to peaks summit")
## dev.off()

