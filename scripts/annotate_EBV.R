
rm(list = ls())

library(GenomicFeatures)
library(GenomicAlignments)
library(GenomicRanges)
library(data.table)
library(org.Hs.eg.db)
library(broom)
library(ggplot2)
library(RColorBrewer)

chrom.info <- read.table(file="/p/keles/SOFTWARE/hg19.chrom.sizes", header=FALSE)
names(chrom.info) <- c('chrom','length')
chrom.info$is_circular <- rep(FALSE, dim(chrom.info)[1])
chrom.info <- chrom.info[-21,]

load("data/ranges/all_EBV_GenomicRanges.RData")

blacklist <- read.table("data/encode/blacklist/hg19-blacklist.txt",header = TRUE)
blacklist <- GRanges(seqnames = blacklist[,1],ranges =
                     IRanges(start = blacklist[,2],end = blacklist[,3]),
                     strand = "*")


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

ranges <- lapply(ranges,set_Seqinfo,chrom.info)

if(!file.exists("inst/gencode/annot.sqlite")){
  txDb <- makeTranscriptDbFromGFF("inst/gencode/gencode.v19.annotation.gtf",format = "gtf",species = "human")
  saveDb(txDb,"inst/gencode/annot.sqlite")
}else{
  txDb <- loadDb("inst/gencode/annot.sqlite")
}

## load genome regions from db
transcript <- transcripts(txDb)
gene <- genes(txDb)
fiveutr <- fiveUTRsByTranscript(txDb)
threeutr <- threeUTRsByTranscript(txDb)
intron <- intronsByTranscript(txDb)
exon <- exonsBy(txDb)

## remove chrM,chrY (since this are female cells) and add chr length info
transcript <- set_Seqinfo(transcript,chrom.info)
gene <- set_Seqinfo(gene,chrom.info)
fiveutr <- set_Seqinfo(unlist(fiveutr),chrom.info)
threeutr <- set_Seqinfo(unlist(threeutr),chrom.info)
intron <- set_Seqinfo(unlist(intron),chrom.info)
exon <- set_Seqinfo(unlist(exon),chrom.info)
blacklist <- set_Seqinfo(unlist(blacklist),chrom.info)

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

## downstream <- function(x,dm=10000)
## {
##   strd <- strand(x)
##   strand(x)[which(strd=='+')] <- '-'
##   strand(x)[which(strd=='-')] <- '+'
##   out <- promoters(x,upstream=dm,downstream=0)
##   strand(out) <- strd
##   return(out)
## }


annoDb <- "org.Hs.eg.db"
gene <- addGeneAnno2feature(gene,annoDb)
## transcript <- addGeneAnno2feature(transcript,annoDb)
## fiveutr <- addGeneAnno2feature(fiveutr,annoDb)
## threeutr <- addGeneAnno2feature(threeutr,annoDb)
## exon <- addGeneAnno2feature(exon,annoDb)
## intron <- addGeneAnno2feature(intron,annoDb)

gene_fwd <- gene[strand(gene) == "+"]
gene_bwd <- gene[strand(gene) == "-"]
chr <- as.character(chrom.info[,1])
gene_fwd <- lapply(chr,function(x,gene)gene[seqnames(gene) == x],gene_fwd)
gene_bwd <- lapply(chr,function(x,gene)gene[seqnames(gene) == x],gene_bwd)
names(gene_fwd) <- chr
names(gene_bwd) <- chr

blacklist <- lapply(chr,function(x,reg)reg[seqnames(reg) == x],blacklist)
names(blacklist) <- chr


promoter_region <- function(gr,dd,stra)
{
  gr <- reduce(gr)
  out <- gr
  if(stra == "+"){
    end(out) <- start(gr) - 1
    start(out) <- end(out) - dd + 1
  }else{
    start(out) <- end(gr) + 1
    end(out) <- start(out) + dd -1
  }
  return(out)
}

prom_fwd <- lapply(gene_fwd,promoter_region,2e3,"+")
prom_bwd <- lapply(gene_bwd,promoter_region,2e3,"-")


overlap_region <- function(prom,gr)
{
  elementMetadata(gr) = NULL
  names(gr) = NULL
  ov <- countOverlaps(prom,gr) > 0
  return(sort(gr[!ov]))
}

gene_fwd <- mapply(overlap_region,prom_fwd,gene_fwd,SIMPLIFY=FALSE)
gene_bwd <- mapply(overlap_region,prom_bwd,gene_bwd,SIMPLIFY=FALSE)

get_probability <- function(glist,chrom.info)
{
  num <- do.call(sum, lapply(glist,width))
  return( num / sum(as.numeric(chrom.info[,2])))
}

prom_prob <- mean(get_probability(prom_fwd,chrom.info),get_probability(prom_bwd,chrom.info))
gene_prob <- mean(get_probability(gene_fwd,chrom.info),get_probability(gene_bwd,chrom.info))

p0 <- c(prom_prob,gene_prob,1 - gene_prob - prom_prob)


region_test <- function(reg,p0)
{
  na_counts <-  sum(is.na(reg$region))
  prom <- sum(reg$region == "5p1",na.rm = TRUE)
  gen <- sum(reg$region == "gene",na.rm = TRUE)
  counts <- c(prom,gen,length(reg) - prom - gen - na_counts)
  p <- counts / sum(counts)
  test <- chisq.test(counts, p = p0)
  labs <- c("Promoter","Gene body","Intergenic")
  dt <- rbind( data.table(lab = labs ,probs = p , case = "ChIP"),
              data.table(lab = labs,probs = p0,case = "genome"))
  plot <- ggplot(dt , aes(case,probs,fill = case))+geom_bar(stat = "identity")+
    scale_fill_brewer(name = "Region",palette = "Pastel1")+facet_grid( . ~ lab)+
    theme(legend.position = "none")+xlab("")+ylab("")
  out <- list( test = test , p = p , data = dt , plot = plot)
  return(out)
}

tests = lapply(ranges,region_test,p0)

whole = lapply(tests,function(x)x$data)
whole = mapply(function(x,y){
  x[,sample:=y]
  return(x)},whole,names(whole),SIMPLIFY=FALSE)
whole = do.call(rbind,whole)
whole = whole[case != "genome"]
whole = rbind(whole , data.table(lab = c("Promoter","Gene body","Intergenic"),
  probs = p0,case = "genome",sample = "genome"))
whole[,sample := factor(sample,levels = c("genome","EBNA2","EBNA3A","EBNA3B","EBNA3C","JK92","JK234","RBPJ"))]

out <-  list(tests,whole)
save(out,  file = "data/RData/annotation.RData" )
## ts <- 5e2
## dp <- 2e3
## um <- 1e4
## dm <- 1e4

## firstexon <- exon[exon$exon_rank==1]
## tss <- promoters(transcript,upstream=ts,downstream=ts)
## promt <- promoters(transcript,upstream=dp,downstream=0)
## upstrm <- promoters(transcript,upstream=um,downstream=0) ##note that upstrm and promt overlaps
## downstrm <- downstream(transcript,dm=dm)

## dis <- 100000
## geneflank <- resize(gene,width=width(gene)+2*dis,fix='center')

## relative_enrichment_test <- function(gr,annot,p = .5,method = "fisher")
## {
##   if(method == "fisher"){
##     ### given two sets of intervals, we need to counts the number of bases in the contigency table
    

    

##   }else if(method == "binom"){
##     ov <- countOverlaps(gr,annot) > 0
##     out <- tidy(binom.test(sum(ov),length(gr),p))
##   }else{
##     stop("Unknown method")
##   }
##   return(out)
## }




## samples <- replicate(1000,sample(length(ranges[[1]]),100),simplify=FALSE)

## p0 = sum(width(promoters(tss,2e3))) / sum(as.numeric(chrom.info[,2]))

## gr <- ranges[[1]]
## annot <- transcript
## p <- mean(countOverlaps(gr,annot) > 0 )

## binom <- mclapply( samples,function(x,gr,annot,p)
##   relative_enrichment_test(gr[x],annot,p,method ="binom"),gr,annot,p, mc.cores = 8)
## binom <- data.table(do.call(rbind,binom))


## ggplot(binom,aes(p.value))+geom_histogram(binwidth = .01, aes(y = ..density..))+
##   geom_vline(xintercept = p,size =2 ,colour = "red",linetype =2)
