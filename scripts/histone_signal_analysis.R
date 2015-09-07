

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

figs_dir <- "figures/enhancers"

p <- ggplot(stats,aes(log10(1 + V1),fill = set))+geom_density(kernel = "rectangular")+
  scale_fill_brewer(palette = "Set1")+facet_grid( histone ~ set , scale = "free_y" )+xlim(0,1.5)+
  theme(legend.position = "none",strip.text.y = element_text(angle = 0),
        axis.text.y = element_text(size = 8),axis.text.x = element_text(size = 8,angle = 45))
                       
pdf(file = file.path(figs_dir,"histone_signal.pdf"))
print(p)
print(p %+% stats_dnase)
dev.off()

## this plots indicates that by filtering into DHS, we are actually getting higher signal, for all the histone marks except H3k27me3, H3k36me3, H3k9me3, H3k4me1 (i.e. the represor marks)

## p1 <- ggplot(stats,aes(histone,V1,color = set))+geom_boxplot()+
##   scale_y_log10()+facet_grid(set ~ .) + scale_color_brewer(palette = "Set1")+
##   theme(legend.position = "none",axis.text.x = element_text(angle = 45))

## pdf(file = file.path(figs_dir,"histone_signal_boxplot.pdf"))
## print(p1)
## print(p1 %+% stats_dnase)
## dev.off()

## says roughly the same as p

## H2A.Z  Histone protein variant (H2A.Z) associated with regulatory elements with dynamic chromatin
## H3K4me1  Mark of regulatory elements associated with enhancers and other distal elements, but also enriched downstream of transcription starts
## H3K4me2  Mark of regulatory elements associated with promoters and enhancers
## H3K4me3  Mark of regulatory elements primarily associated with promoters/transcription starts
## H3K9ac  Mark of active regulatory elements with preference for promoters
## H3K9me1  Preference for the 5′ end of genes
## H3K9me3  Repressive mark associated with constitutive heterochromatin and repetitive elements
## H3K27ac  Mark of active regulatory elements; may distinguish active enhancers and promoters from their inactive counterparts
## H3K27me3 Repressive mark established by polycomb complex activity associated with repressive domains and silent developmental genes
## H3K36me3 Elongation mark associated with transcribed portions of genes, with preference for 3′ regions after intron 1
## H3K79me2  Transcription-associated mark, with preference for 5′ end of genes
## H4K20me1  Preference for 5′ end of genes

## following table 2 from Nature's first encode consortia paper
## http://www.nature.com/nature/journal/v489/n7414/full/nature11247.html
## (which is above) we can remove from our list:

## H2A.Z, H3K79me2 , H4K20me1 , the rest are kept to coincide with fig 3C


histones <- unique(stats[,(histone)])
histones <- grepV("H3",histones)
histones <- grepV("79",histones,invert = TRUE)


create_matrix <- function(set,data,cols)
{
  setkey(data,histone)
  data <- copy(data[cols])
  setkey(data , set)
  data <- data[set]
  data <- dcast.data.table(data =data , match ~ histone ,value.var = "V1")
  data <- log10(1 + as.matrix(data[,-1,with = FALSE]))
  return(data)
}


heatmap_analysis <- function(set, data , cols)
{
  data <- create_matrix(set,data,cols)
  heatmap(data , col = brewer.pal(11,"RdYlGn"),labRow = "",main = set)
}


pdf(file = file.path(figs_dir,"EBNA2_heatmaps.pdf"),height = 8)
heatmap_analysis("EBNA2",stats,histones)
heatmap_analysis("EBNA2",stats_dnase,histones)
dev.off()

pdf(file = file.path(figs_dir,"EBNA3A_heatmaps.pdf"),height = 8)
heatmap_analysis("EBNA3A",stats,histones)
heatmap_analysis("EBNA3A",stats_dnase,histones)
dev.off()

pdf(file = file.path(figs_dir,"EBNA3B_heatmaps.pdf"),height = 8)
heatmap_analysis("EBNA3B",stats,histones)
heatmap_analysis("EBNA3B",stats_dnase,histones)
dev.off()

pdf(file = file.path(figs_dir,"EBNA3C_heatmaps.pdf"),height = 8)
heatmap_analysis("EBNA3C",stats,histones)
heatmap_analysis("EBNA3C",stats_dnase,histones)
dev.off()

pdf(file = file.path(figs_dir,"RBPJ_heatmaps.pdf"),height = 8)
heatmap_analysis("RBPJ",stats,histones)
heatmap_analysis("RBPJ",stats_dnase,histones)
dev.off()

## the clustering pattern for the histone marks is topologically equivalente among
## all datasets. It may be worth to look into the LDA idea of normal clustering

pam_analysis <- function(set,K,data,cols)
{
  data <- create_matrix(set,data,cols)
  pam_analysis <- pam(data,K)

  out <- list()
  out[["nr_clusters"]] <- table(pam_analysis$clustering)

  medoids <- data.table(pam_analysis$medoids)
  medoids[,cluster := 1:K]

  out[["medoids"]] <- medoids  
  medoids <- melt(medoids,id.vars = "cluster")


  rf <- colorRampPalette(brewer.pal(11,"RdYlGn"))

  out[["plot"]] <- ggplot(medoids , aes(variable , as.factor(cluster) ,fill = value))+
    geom_tile()+ xlab("Histone marks")+
    scale_fill_gradientn(colours = rf(8))+
    theme(legend.position = "top",axis.text.x = element_text(angle = 90))+
    ggtitle(set)

  return(out)
}

K <- 10

result <- list()
result_dnase <- list()

for(set in sets){
  pdf(file = file.path(figs_dir,paste0(set,"_medoids_K",K,".pdf")))
  result[[set]] <- pam_analysis(set,K,stats,histones)
  result_dnase[[set]] <- pam_analysis(set,K,stats_dnase,histones)
  print(result[[set]][["nr_clusters"]])
  print(result[[set]][["plot"]])
  print(result_dnase[[set]][["nr_clusters"]])
  print(result_dnase[[set]][["plot"]])
  dev.off()
}


results <- list()
results[["all"]] <- result
results[["dnase"]] <- result_dnase

save(results,file = "data/RData/pam_analysis_K10.RData")
