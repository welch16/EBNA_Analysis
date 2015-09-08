

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

set_dnase <- lapply(sets,
  function(set,peak_list)peak_list[[set]][Dnase == 1,(match)],
    peak_list)
names(set_dnase) <- sets

profiles_dnase <- lapply(sets,
  function(ss,set_dnase,stats){
    matchs <- set_dnase[[ss]]
    out <- stats[[ss]]
    out <- lapply(out,function(x,matchs){
      setkey(x,match)
      return(x[matchs])
    },matchs)
    return(out)},set_dnase,profiles)
names(profiles_dnase) <- sets

trim <- .05
suppressWarnings(
histone_profiles <- lapply(sets,
  function(set,profiles,t){
    out <- profiles[[set]]
    out <- mclapply(out,function(x,t)x[,mean(tagCounts,trim = t, na.rm = TRUE),by = coord],
        t,mc.cores = mc)
    return(out)},profiles,trim)
)
names(histone_profiles) <- sets
           
suppressWarnings(
histone_profiles_dnase <- lapply(sets,
  function(set,profiles,t){
    out <- profiles[[set]]
    out <- mclapply(out,function(x,t)x[,mean(tagCounts,trim = t, na.rm = TRUE),by = coord],
        t,mc.cores = mc)
    return(out)},profiles_dnase,trim)
)                    
       
names(histone_profiles_dnase) <- sets

histone_profiles <- lapply(sets,add_names,histone_profiles)
histone_profiles <- lapply(histone_profiles,function(x)x[rep == 1])
histone_profiles <- do.call(rbind,histone_profiles)
histone_profiles[,rep := NULL]

histone_profiles_dnase <- lapply(sets,add_names,histone_profiles_dnase)
histone_profiles_dnase <- lapply(histone_profiles_dnase,function(x)x[rep == 1])
histone_profiles_dnase <- do.call(rbind,histone_profiles_dnase)
histone_profiles_dnase[,rep := NULL]

histones <- histone_profiles[,unique(histone)]
histones <- grepV("H3",histones)
histones <- grepV("79",histones,invert = TRUE)

setkey(histone_profiles,histone)
setkey(histone_profiles_dnase,histone)
histone_profiles <- histone_profiles[histones]
histone_profiles_dnase <- histone_profiles_dnase[histones]

figs_dir <- "figures/for_paper"

setkey(histone_profiles,set)
setkey(histone_profiles_dnase,set)

pdf(file = file.path(figs_dir, "fig3C_EBNA3B.pdf")  ,width = 9,height = 8)
ggplot(histone_profiles["EBNA3B"],aes(coord,V1,colour = histone))+geom_line(size = 1.2)+
  scale_color_brewer(palette = "Dark2")+
  theme(legend.position = "bottom")+geom_vline(xintercept = 0 , linetype =2 ,size = 1.2)+
  guides(colour = guide_legend(nrow = 2))+ylim(0,2)+ylab("Normalized signal")+
  xlab("Distance to peaks summit")
ggplot(histone_profiles_dnase["EBNA3B"],aes(coord,V1,colour = histone))+geom_line(size = 1.2)+
  scale_color_brewer(palette = "Dark2")+
  theme(legend.position = "bottom")+geom_vline(xintercept = 0 , linetype =2 ,size = 1.2)+
  guides(colour = guide_legend(nrow = 2))+ylim(0,2)+ylab("Normalized signal")+
  xlab("Distance to peaks summit")
dev.off()

pdf(file = file.path(figs_dir, "fig3C_EBNA3A.pdf")  ,width = 9,height = 8)
ggplot(histone_profiles["EBNA3A"],aes(coord,V1,colour = histone))+geom_line(size = 1.2)+
  scale_color_brewer(palette = "Dark2")+
  theme(legend.position = "bottom")+geom_vline(xintercept = 0 , linetype =2 ,size = 1.2)+
  guides(colour = guide_legend(nrow = 2))+ylim(0,2)+ylab("Normalized signal")+
  xlab("Distance to peaks summit")
ggplot(histone_profiles_dnase["EBNA3A"],aes(coord,V1,colour = histone))+geom_line(size = 1.2)+
  scale_color_brewer(palette = "Dark2")+
  theme(legend.position = "bottom")+geom_vline(xintercept = 0 , linetype =2 ,size = 1.2)+
  guides(colour = guide_legend(nrow = 2))+ylim(0,2)+ylab("Normalized signal")+
  xlab("Distance to peaks summit")
dev.off()

pdf(file = file.path(figs_dir, "fig3C_EBNA3C.pdf")  ,width = 9,height = 8)
ggplot(histone_profiles["EBNA3C"],aes(coord,V1,colour = histone))+geom_line(size = 1.2)+
  scale_color_brewer(palette = "Dark2")+
  theme(legend.position = "bottom")+geom_vline(xintercept = 0 , linetype =2 ,size = 1.2)+
  guides(colour = guide_legend(nrow = 2))+ylim(0,2)+ylab("Normalized signal")+
  xlab("Distance to peaks summit")
ggplot(histone_profiles_dnase["EBNA3C"],aes(coord,V1,colour = histone))+geom_line(size = 1.2)+
  scale_color_brewer(palette = "Dark2")+
  theme(legend.position = "bottom")+geom_vline(xintercept = 0 , linetype =2 ,size = 1.2)+
  guides(colour = guide_legend(nrow = 2))+ylim(0,2)+ylab("Normalized signal")+
  xlab("Distance to peaks summit")
dev.off()

pdf(file = file.path(figs_dir, "fig3C_EBNA2.pdf")  ,width = 9,height = 8)
ggplot(histone_profiles["EBNA2"],aes(coord,V1,colour = histone))+geom_line(size = 1.2)+
  scale_color_brewer(palette = "Dark2")+
  theme(legend.position = "bottom")+geom_vline(xintercept = 0 , linetype =2 ,size = 1.2)+
  guides(colour = guide_legend(nrow = 2))+ylim(0,2)+ylab("Normalized signal")+
  xlab("Distance to peaks summit")
ggplot(histone_profiles_dnase["EBNA2"],aes(coord,V1,colour = histone))+geom_line(size = 1.2)+
  scale_color_brewer(palette = "Dark2")+
  theme(legend.position = "bottom")+geom_vline(xintercept = 0 , linetype =2 ,size = 1.2)+
  guides(colour = guide_legend(nrow = 2))+ylim(0,2)+ylab("Normalized signal")+
  xlab("Distance to peaks summit")
dev.off()

pdf(file = file.path(figs_dir, "fig3C_RBPJ.pdf")  ,width = 9,height = 8)
ggplot(histone_profiles["RBPJ"],aes(coord,V1,colour = histone))+geom_line(size = 1.2)+
  scale_color_brewer(palette = "Dark2")+
  theme(legend.position = "bottom")+geom_vline(xintercept = 0 , linetype =2 ,size = 1.2)+
  guides(colour = guide_legend(nrow = 2))+ylim(0,2)+ylab("Normalized signal")+
  xlab("Distance to peaks summit")
ggplot(histone_profiles_dnase["RBPJ"],aes(coord,V1,colour = histone))+geom_line(size = 1.2)+
  scale_color_brewer(palette = "Dark2")+
  theme(legend.position = "bottom")+geom_vline(xintercept = 0 , linetype =2 ,size = 1.2)+
  guides(colour = guide_legend(nrow = 2))+ylim(0,2)+ylab("Normalized signal")+
  xlab("Distance to peaks summit")
dev.off()

