

rm(list = ls())
library(GenomicAlignments)
library(ggplot2)
library(data.table)
library(ChIPUtils)
library(GenomicFeatures)
library(org.Hs.eg.db)


## load West peaks

west_dr <- "./data/review_comments"
west_peaks <- list.files(west_dr)

west <- lapply(file.path(west_dr,west_peaks),read.table,skip = 1)
west <- lapply(west,data.table)
names(west) <- c("EBNA2","EBNA3C")

load(file = "data/RData/unified_lists_wProbs.RData") ## unified_lists
load("data/RData/factor_overlaps.RData")

DT <- cbind(unified_lists$peaks, unified_lists$overlaps)

DT <- cbind(DT,factor_overlaps)

peaks <- dt2gr(DT[,1:3,with = FALSE])

west <- lapply(west,function(x){
  setnames(x,names(x),c("seqnames","start","end","name","score"))
  return(x)})

west_ebna2 <- ifelse(countOverlaps(peaks,dt2gr(west[["EBNA2"]])) > 0,1,0)
west_ebna3C <- ifelse(countOverlaps(peaks,dt2gr(west[["EBNA3C"]])) > 0,1,0)

DT[,"West.EBNA2" := west_ebna2]
DT[,"West.EBNA3C" := west_ebna3C]



if(!file.exists("inst/gencode/annot.sqlite")){
  txDb <- makeTranscriptDbFromGFF("inst/gencode/gencode.v19.annotation.gtf",
                                  format = "gtf",species = "human")
  saveDb(txDb,"inst/gencode/annot.sqlite")
}else{
  txDb <- loadDb("inst/gencode/annot.sqlite")
}

gene <- genes(txDb)

## remove chrM,chrY (since this are female cells) and add chr length info

chrom.info <- read.table(file="/p/keles/SOFTWARE/hg19.chrom.sizes", header=FALSE)
names(chrom.info) <- c('chrom','length')
chrom.info$is_circular <- rep(FALSE, dim(chrom.info)[1])
chrom.info <- chrom.info[-21,]

set_Seqinfo <- function(gr,chrom.info)
{
  chr <- as.character(chrom.info[,1])
  seqlevels(gr,force = TRUE) <- chr
  idx <- match(seqnames(seqinfo(gr)),chr)
  ll <- chrom.info[,2]
  names(ll) <- chr
  yy <- chrom.info[,3]
  names(yy) <-  chr
  seqlengths(gr) <- ll[idx]
  seqinfo(gr)@is_circular <- yy[idx]
  return(gr)
}

gene <- set_Seqinfo(gene,chrom.info)


## blacklist <- set_Seqinfo(unlist(blacklist),chrom.info)


getFirstHitIndex <- function(x)sapply(unique(x), function(i) which(x == i)[1])


addGeneAnno <- function(annoDb, geneID, type= "Entrez Gene ID")
{
  kk <- unlist(geneID)
  require(annoDb, character.only = TRUE)
  annoDb <- eval(parse(text=annoDb))
  if (type == "Entrez Gene ID") {
    kt <- "ENTREZID"
  } else if (type =="Ensembl gene ID" || type == "Ensembl Gene ID") {
    kt <- "ENSEMBL"
  } else {
    warnings("geneID type is not supported...\tPlease report it to developer...\n")
    return(NA)
  }
  if (sum(kk %in% keys(annoDb, "ENSEMBL")) > 0){
    ann <- suppressWarnings(select(annoDb,keys=kk,
                keytype=kt,
                columns=c("ENTREZID", "ENSEMBL", "SYMBOL", "GENENAME")))
    idx <- getFirstHitIndex(ann[,kt])
    ann <- ann[idx,]
    idx <- unlist(sapply(kk, function(x) which(x==ann[,kt])))
    ann <- ann[idx,]
  } else {
    ann <- matrix(NA, ncol = 4, nrow = length(geneID))
    colnames(ann) <- c("ENTREZID", "ENSEMBL", "SYMBOL", "GENENAME")
  }
  return(ann)
}

addGeneAnno2feature <- function(feature,annoDb,type='Ensembl Gene ID')
{
  gid <- sapply(feature$gene_id, function(x) strsplit(x, split='.',fix=TRUE)[[1]][1])
  geneAnno <- addGeneAnno(annoDb, gid,type=type)
  for(ii in names(geneAnno)){
    elementMetadata(feature)[[ii]] <- geneAnno[[ii]]
  }
  return(feature)
}

annoDb <- "org.Hs.eg.db"
gene <- addGeneAnno2feature(gene,annoDb)

gene_fwd <- gene[strand(gene) == "+"]
gene_bwd <- gene[strand(gene) == "-"]
chr <- as.character(chrom.info[,1])

gene_fwd <- lapply(chr,function(x,gene)gene[seqnames(gene) == x],gene_fwd)
gene_bwd <- lapply(chr,function(x,gene)gene[seqnames(gene) == x],gene_bwd)

names(gene_fwd) <- chr
names(gene_bwd) <- chr

aux_peaks <- as.list(split(peaks,seqnames(peaks)))[chr]

closest_distance_gene <- function(peak,fwd,bwd)
{
  fwd_dist <- distanceToNearest(peak,fwd,ignore.strand = TRUE)
  bwd_dist <- distanceToNearest(peak,bwd,ignore.strand = TRUE)

  fwd_symbol <- mcols(fwd)[subjectHits(fwd_dist),"SYMBOL"]
  bwd_symbol <- mcols(bwd)[subjectHits(bwd_dist),"SYMBOL"]

  fwd_name <- mcols(fwd)[subjectHits(fwd_dist),"GENENAME"]
  bwd_name <- mcols(bwd)[subjectHits(bwd_dist),"GENENAME"]

  fwd_ens <- mcols(fwd)[subjectHits(fwd_dist),"ENSEMBL"]
  bwd_ens <- mcols(bwd)[subjectHits(bwd_dist),"ENSEMBL"]

  
  fdist <- mcols(fwd_dist)$distance
  bdist <- mcols(bwd_dist)$distance

  fwd <- data.table(fwd_nearest_dist = fdist, fwd_symbol = fwd_symbol,
                    fwd_name = fwd_name,fwd_ensembl = fwd_ens)
  bwd <- data.table(bwd_nearest_dist = bdist, bwd_symbol = bwd_symbol,
                    bwd_name = bwd_name, bwd_ensembl = bwd_ens)

  dd = gr2dt(peak)
  dd[, match := paste0(seqnames,":",start,"-",end)]

  out <- cbind(fwd,bwd)

  out[,match := dd$match]
  
  return(out)  
}

annot_genes <- mapply(closest_distance_gene,aux_peaks,gene_fwd,gene_bwd,
  SIMPLIFY = FALSE)

annot_genes <- do.call(rbind,annot_genes)

DT[, match := paste0(seqnames,":",start,"-",end)]                      

DT <- merge(DT,annot_genes,by = "match")

DT[,match :=  NULL]
write.table(DT ,"inst/west/EBV_peaks_West_ov_gene_annot.csv",
   sep = "\t",quote = FALSE,col.names = TRUE,row.names = FALSE)           


