

rm(list = ls())
library(GenomicAlignments)
library(ggplot2)
library(data.table)
library(ChIPUtils)

## load West peaks

west_dr <- "./data/review_comments"
west_peaks <- list.files(west_dr)

west <- lapply(file.path(west_dr,west_peaks),read.table,skip = 1)
west <- lapply(west,data.table)
names(west) <- c("EBNA2","EBNA3C")

load("data/ranges/all_EBV_GenomicRanges.RData")

west <- lapply(west,function(x){
  setnames(x,names(x),c("seqnames","start","end","name","score"))
  return(x)})

west_gr <- lapply(west,function(x)dt2gr(x[,1:3,with =FALSE]))

ranges <- lapply(ranges,function(x){
  mcols(x)$West_EBNA2 <- ifelse(
    countOverlaps(x,west_gr[["EBNA2"]]) > 0,1,0)
  return(x)})

ranges <- lapply(ranges,function(x){
  mcols(x)$West_EBNA3C <- ifelse(
    countOverlaps(x,west_gr[["EBNA3C"]]) > 0,1,0)
  return(x)})

ranges_dt <- lapply(ranges,gr2dt)

west_analysis <- function(dt,nm,which,what,dnase = TRUE)
{                          
  stopifnot(which %in% c("EBNA2","EBNA3C"))
  stopifnot(what %in% c("tf","histone","all","ours"))

  nms <- names(dt)
  nms <- nms[-c(1:8)]
  nms <- nms[ nms != "allTF"]

  west <- nms[grep("West",nms)]
  nms <- nms[grep("West",nms,invert = TRUE)]

  ours <- c("EBNA2","EBNA3A","EBNA3B","EBNA3C","JK92","JK234","RBPJ")
  nms <- nms[-sapply(ours,grep,nms)]
  if(dnase){
    dt <- dt[Dnase == 1]
  }
  nms <- nms[grep("Dnase",nms,invert = TRUE)]
  
  histone <- c("H2ZA","H3K27AC","H3K27ME3","H3K36ME3","H3K4ME1",
               "H3K4ME2","H3K4ME3","H3K79ME2","H3K9AC","H3K9ME3",
               "H4K20ME1")
  nms <- nms[-sapply(histone,grep,nms)]
  if(what == "all"){
    cols <- c("ours","histone","tf")
  }else if(what == "ours"){
    cols <- ours
  }else if(what == "tf"){
    cols <- nms
  }else{
    cols <- histone
  }
  which <- west[grep(which,west)]
  idx <- dt[[which]] == 1

    browser()


  dt <- dt[,cols,with = FALSE]

  props <- split(dt,idx)
  props <- lapply(props,colMeans)
  props <- lapply(props,function(x){
    nn <- names(x)
    names(x) <- NULL
    data.table(col = nn,prop = x)})
  props <- mapply(function(x,y)x[,overlap := y],props,names(props),
                  SIMPLIFY = FALSE)
  props <- do.call(rbind,props)
  
  ## get indx
  ggplot(props,aes(col,prop,fill = overlap))+
    geom_bar(stat ="identity",position = "dodge",width = .5)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90))+
    scale_fill_brewer(palette = "Set1")+ggtitle(nm)

}



pdf(width = 14,height = 6)
mapply(west_analysis,ranges_dt,names(ranges_dt),
       MoreArgs = list("EBNA3C","tf"),SIMPLIFY = FALSE)
dev.off()






## colors = c("blue","red","chartreuse4","darkorchid4","darkorange","goldenrod4","yellow")
## names(colors) = c("EBNA3A","EBNA3B","EBNA3C","EBNA2","RBPJ","JK234","JK92")

## # Get EBNA3 range
## rangeDir <- "../Generated/RData"
## figsDir <-  "../Generated/Figs"
## load(file = file.path(rangeDir,"all_EBV_GenomicRanges.RData"))  # ranges
## aux = c(ranges[["EBNA3A"]],ranges[["EBNA3B"]],ranges[["EBNA3C"]])
## ebna3 = reduce(GRanges(seqnames = seqnames(aux),ranges= IRanges(start = start(aux),end = end(aux)),strand = "*"))





## # Get West peaks
## westDir <-  "../WestMJ_peaks"
## file = list.files(westDir)

## west_peaks = read.table(file.path(westDir,file),skip = 1,
##   colClasses = c("character","numeric","numeric","character","numeric"))

## west_range = GRanges(seqnames = west_peaks[[1]],ranges = IRanges(start = west_peaks[[2]],
##   end = west_peaks[[3]]),strand = "*")

## mg = c(0,100,200,500,1000,2000)

## overlap.weights <-function(myRange,range,maxGap)
## {
##   n = length(myRange)
##   m = length(range)
##   interCount = sum(countOverlaps(myRange,range,maxgap = maxGap))
##   return(c("11"=interCount,"10"=n-interCount,"01"=m-interCount))
## }

## # For EBNA3A
## ebna3a_venn = lapply(mg,function(x,ranges,west_range)Venn(SetNames = c(paste0("our_ebna3a_mg=",x),"west"),
##   Weight = overlap.weights(ranges[["EBNA3A"]],west_range,x)),
##   ranges,west_range)
## names(ebna3a_venn) =mg

## pdf(file =file.path(figsDir,"ebna3a_vs_west_venn.pdf"))
## lapply(ebna3a_venn,FUN = plot)
## dev.off()

## # For EBNA3B
## ebna3b_venn = lapply(mg,function(x,ranges,west_range)Venn(SetNames = c(paste0("our_ebna3b_mg=",x),"west"),
##   Weight = overlap.weights(ranges[["EBNA3B"]],west_range,x)),
##   ranges,west_range)
## names(ebna3b_venn) =mg

## pdf(file =file.path(figsDir,"ebna3b_vs_west_venn.pdf"))
## lapply(ebna3b_venn,FUN = plot)
## dev.off()

## # For EBNA3C
## ebna3c_venn = lapply(mg,function(x,ranges,west_range)Venn(SetNames = c(paste0("our_ebna3c_mg=",x),"west"),
##   Weight = overlap.weights(ranges[["EBNA3C"]],west_range,x)),
##   ranges,west_range)
## names(ebna3c_venn) =mg

## pdf(file =file.path(figsDir,"ebna3c_vs_west_venn.pdf"))
## lapply(ebna3c_venn,FUN = plot)
## dev.off()


## # For EBNA3C
## ebna3_venn = lapply(mg,function(x,ranges,west_range)Venn(SetNames = c(paste0("our_ebna3_mg=",x),"west"),
##   Weight = overlap.weights(ranges,west_range,x)),
##   ebna3,west_range)
## names(ebna3_venn) =mg

## pdf(file =file.path(figsDir,"ebna3_all_vs_west_venn.pdf"))
## lapply(ebna3_venn,FUN = plot)
## dev.off()

