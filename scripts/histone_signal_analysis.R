

rm( list = ls())

library(data.table)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(parallel)
library(pracma)

profile_dir <- "data/RData/profile_matrices"

sets <- c("EBNA2","EBNA3A","EBNA3B","EBNA3C","RBPJ")

files <- list.files(profile_dir)
window_ext <- 2000
mc <- detectCores()

load("data/RData/unified_lists_wProbs.RData")
load("data/RData/factor_overlaps.RData")

load_clean_profile <- function(file ,window_ext)
{
  load(file) ## profile
  profile[,match := paste(chr,match,sep = "_")]
  profile[,chr := NULL]
  profile <- profile[order(match,coord)]
  npeaks <- nrow(profile) / (2 * window_ext + 1)
  profile[,coord := rep(-window_ext:window_ext, npeaks)]
  return(profile)
}


grepV <- function(pattern,x,ignore.case = FALSE,perl = FALSE,value = FALSE,
  fixed = FALSE,useBytes = FALSE, invert = FALSE)
{
  return(x[grep(pattern,x,ignore.case = FALSE,perl = FALSE,value = FALSE,
  fixed = FALSE,useBytes = FALSE, invert = FALSE)])
}


files1 <- grepV("EBNA3A",files)
files2 <- grepV("EBNA3B",files)
files3 <- grepV("EBNA3C",files)
files4 <- grepV("EBNA2",files)
files5 <- grepV("RBPJ",files)

profiles <- list()
profiles[["EBNA3A"]] <- mclapply(file.path(profile_dir,files1),
  load_clean_profile,window_ext, mc.cores = mc,mc.preschedule = TRUE)

profiles[["EBNA3B"]] <- mclapply(file.path(profile_dir,files2),
  load_clean_profile,window_ext, mc.cores = mc,mc.preschedule = TRUE)

profiles[["EBNA3C"]] <- mclapply(file.path(profile_dir,files3),
  load_clean_profile,window_ext, mc.cores = mc,mc.preschedule = TRUE)

profiles[["EBNA2"]] <- mclapply(file.path(profile_dir,files4),
  load_clean_profile,window_ext, mc.cores = mc,mc.preschedule = TRUE)

profiles[["RBPJ"]] <- mclapply(file.path(profile_dir,files5),
  load_clean_profile,window_ext, mc.cores = mc,mc.preschedule = TRUE)


add_names <- function(set,data,files,mc = 24)
{
    out <- mcmapply(function(x,y){
      v <- strsplit(y,"_")[[1]]
      aux <- copy(x)
      aux[,set := v[1]]
      aux[,histone := v[2]]
      aux[,rep := as.numeric(gsub("rep","",v[3]))]
      return(aux)
    },data[[set]],grepV(set,files),SIMPLIFY = FALSE,mc.cores = mc)
  
  return(out)
}

test <- list()
test[["RBPJ"]] <- mclapply(profiles[["RBPJ"]],
  function(x){    
    out <- x[,mean(tagCounts,trim = .05,na.rm = TRUE),by = coord]
    return(out)},mc.cores = mc)             
test <- add_names("RBPJ",test,files)
test <- do.call(rbind,test)


pdf(file = "RBPJ.pdf")
ggplot(test[grep("H3",histone)],
       aes(coord,V1,colour = histone,linetype = as.factor(rep)))+
  geom_line(linesize = 1)+scale_color_brewer(palette = "Set1")+
  theme(panel.background = element_rect(fill = "black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ylim(0,2)
dev.off()


M <- 1.5e3
lb <- -M
ub <- M
stats <- lapply(profiles,function(x,lb,ub){
  out <- mclapply(x,function(y,lb,ub){
    y[between(coord,lb,ub),max(tagCounts)/min(tagCounts),by = match]
  },lb,ub,mc.cores = mc)
  return(out)},lb,ub)


stats <- lapply(stats,function(x)do.call(rbind,x))
stats <- do.call(rbind,stats)

ggplot(stats[,mean(signal),by = .(match,set,histone)] ,
       aes(histone , V1))+geom_boxplot()+
  facet_grid( set ~ .)+scale_y_log10()+ylab("signal")+xlab("histone")+
  theme(axis.text.x = element_text(angle = 90))
dev.off()
