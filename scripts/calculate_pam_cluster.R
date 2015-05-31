
rm(list = ls())
graphics.off()

library(data.table)
library(cluster)
library(GenomicRanges)




## loads binary matrices + functions to generate them
load("data/ranges/all_EBV_GenomicRanges.RData")

source("R/heatmap_functions.R")
source("R/clustering_analysis.R")


peaks = c("EBNA2","EBNA3A","EBNA3B","EBNA3C","RBPJ")
pattern = c("EBNA","RBPJ","H?K","Dnase","H2Z")

cols_to_remove <- function(patterns,set)
{
  ll = lapply(pattern,FUN = grepl_cmd,set)
  n = length(ll)
  ll = as.list(paste0(ll,c(rep("&",n-1),"")))
  return(do.call(paste0,ll))
}

get_names <- function(range_data) return(names(range_data@elementMetadata@listData))

grepl_cmd <- function(variable,set)
{
return(paste0("!grepl(",change_to_str_val(variable),",",set,")",sep =""))
}

change_to_str_val <- function(str) return(paste0("'",str,"'"))

mc = 8


mats <- mcmapply(build_binary_matrix, names(ranges),ranges,
  MoreArgs = list(expr = "Dnase == 1",
             col_expr = cols_to_remove(pattern,'columns')),
			 SIMPLIFY=FALSE,mc.silent=TRUE,mc.cores = mc)

maxClusters <- 15
clusters <- lapply(mats[-(5:6)],function(x,maxClusters){
  mclapply(1:(maxClusters-1),function(k,y){
    pam_silhouette_dt(y,k+1)},x,mc.cores = mc,mc.preschedule = TRUE,
           mc.silent = TRUE)},maxClusters)

save(clusters,file = "data/RData/pam_clusters.RData")
