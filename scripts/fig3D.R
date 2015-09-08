
rm( list = ls())

library(data.table)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(parallel)
library(pracma)
library(cluster)
library(scales)

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


M <- 250
lb <- -M
ub <- M

stats <- lapply(profiles,function(x,lb,ub){
  out <- mclapply(x,function(y,lb,ub){
    y[between(coord,lb,ub),max(tagCounts),by  = match]
  },lb,ub,mc.cores = mc)
  return(out)},lb,ub)

## stats <- lapply(profiles,function(x,lb,ub){
##   out <- mclapply(x,function(y,lb,ub){
##     y[between(coord,lb,ub),trapz(coord,tagCounts),by  = match]
##   },lb,ub,mc.cores = mc)
##   return(out)},lb,ub)

stats <- lapply(sets,add_names,stats)
names(stats) <- sets
stats <- do.call(rbind,stats)

stats <- stats[rep == 1]
stats[,rep := NULL]

set_dnase <- lapply(sets,
  function(set,peak_list)peak_list[[set]][Dnase == 1,(match)],
    peak_list)
names(set_dnase) <- sets

stats_dnase <- lapply(sets,
  function(ss,set_dnase,stats){
    out <- stats[set == ss]
    matchs <- set_dnase[[ss]]
    setkey(out,match)
    return(out[matchs])},set_dnase,stats)
stats_dnase <- do.call(rbind,stats_dnase)

setkey(stats,histone)
setkey(stats_dnase,histone)

signal <- stats["H3k27ac"]
signal_dnase <- stats_dnase["H3k27ac"]

signal <- signal[order(set,V1)]
signal_dnase <- signal_dnase[order(set,V1)]

lengths <- signal[,length(V1),by = set]
lengths_dnase <- signal_dnase[,length(V1), by = set]

signal[, rank := lengths[,1:V1, by = set][,(V1)] ]
signal_dnase[,rank := lengths_dnase[,1:V1, by =set][,(V1)] ]

setkey(signal,set)
setkey(signal_dnase,set)

slopes <- signal[ , (max(V1) - min(V1))/(max(rank) - min(rank)) , by = set]
slopes_dnase <- signal_dnase[ , (max(V1) - min(V1))/(max(rank) - min(rank)) , by = set]

setnames(slopes,names(slopes),c("set","slope"))
setnames(slopes_dnase,names(slopes_dnase),c("set","slope"))

get_rank <- function(set,signal,slopes,fold)
{
  diff <- diff(signal[set,(V1)])
  slope <- slopes[set, (slope)]
  mean <- mean(signal[set , (V1)])
  candidates <- which(abs(diff - slope ) < 1e-1)
  min_rank <- min(signal[set][candidates][V1  / mean > fold][,(rank)])
  return(min_rank)
}

build_rect <- function(set,signal,slopes,fold)
{
  rank <- get_rank(set,signal,slopes,fold)
  line <- signal[set][rank]
  low <- data.table(xmin = 0,xmax = line[,(rank)],ymin = 0,ymax = line[,(V1)],fill = "z")
  high <- data.table(xmin = line[,(rank)],xmax = Inf , ymin = line[,(V1)],ymax = Inf, fill ="a")
  lims <- rbind(low,high)
  return(lims)
}

figs_dir <- "figures/for_paper"

fold <- 4

for(set in sets){
  lims <- build_rect(set,signal,slopes,fold)
  lims_dnase <- build_rect(set,signal_dnase,slopes_dnase,fold)
  pdf(file = file.path(figs_dir,paste0("fig3D_",set,".pdf")))
  p <- ggplot(signal[set] , aes(rank , V1 ,colour = set))+geom_line()+
    scale_colour_brewer(palette = "Set1")+
    theme(legend.position = "none")+ylab("H3K27ac signal")+
    geom_rect(data = lims,aes(xmin = xmin,xmax = xmax,ymin = ymin ,ymax = ymax, fill = fill),
      alpha = .2,inherit.aes = FALSE)+scale_fill_brewer(palette = "Set1")+
    geom_abline(slope = 0,intercept = min(lims[,(ymax)]),linetype =2,size = .2)+
    geom_vline(xintercept = min(lims[,(xmax)]),linetype = 2,size = .2)
  print(p)
  p1 <- ggplot(signal_dnase[set] , aes(rank , V1,colour = set))+geom_line()+
    scale_colour_brewer(palette = "Set1")+
    theme(legend.position = "none")+ylab("H3K27ac signal")+
      geom_rect(data = lims_dnase,aes(xmin = xmin,xmax = xmax,ymin = ymin ,ymax = ymax, fill = fill),
       alpha = .2, inherit.aes = FALSE)+scale_fill_brewer(palette = "Set1")+
    geom_abline(slope = 0,intercept = min(lims_dnase[,(ymax)]),linetype =2,size = .2)+
    geom_vline(xintercept = min(lims_dnase[,(xmax)]),linetype = 2,size = .2)
  print(p1)
  dev.off()
}
