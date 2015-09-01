
rm(list = ls())

library(data.table)
library(VennDiagram)
library(GenomicRanges)

load("data/RData/unified_lists_wProbs.RData")
load("data/RData/factor_overlaps.RData")

west <- data.table(read.table("inst/west/WestMJ_GSM1153765_EBNA3C_peaks.txt",
                              header = FALSE,skip = 1))
setnames(west, names(west),c("seqnames","start","end","name","score"))
west[,name := NULL]

peaks <- cbind(unified_lists$peaks,unified_lists$overlaps,Dnase = factor_overlaps$Dnase)
peaks[,width := NULL]


dt2gr <- function(dt)
{
  GRanges(seqnames = dt[,(seqnames)],
          ranges = IRanges(start = dt[,(start)],
            end = dt[,(end)]),strand = "*")
}

mg <- c(0,100,200,500,1000,2000)

venn_analysis <- function( peaks , west ,mg,name)
{
  ov <- findOverlaps( dt2gr(peaks),dt2gr(west),maxgap = mg)
  out <- draw.pairwise.venn(area1 = queryLength(ov),
                            area2 = subjectLength(ov),
                            cross.area = length(ov),
                            category = c("West",name),
                            cex = 1.8,
                            lwd = 3,
                            cat.cex = 1.8,
                            col = c("red","blue"),
                            scale = FALSE)
  return(out)
}



figs_dir <- "figures/for_paper"

pdf(file = file.path(figs_dir,"figS1_3A.pdf"))
venn <- venn_analysis(peaks[ EBNA3A == 1],west, 0, "EBNA3A")
grid.draw(venn)
grid.newpage()
venn <- venn_analysis(peaks[ EBNA3A == 1 & Dnase == 1],west, 0, "EBNA3A\nand\nDHS")
grid.draw(venn)    
dev.off()


pdf(file = file.path(figs_dir,"figS1_3B.pdf"))
venn <- venn_analysis(peaks[ EBNA3B == 1],west, 0, "EBNA3B")
grid.draw(venn)
grid.newpage()
venn <- venn_analysis(peaks[ EBNA3B == 1 & Dnase == 1],west, 0, "EBNA3B\nand\nDHS")
grid.draw(venn)    
dev.off()


pdf(file = file.path(figs_dir,"figS1_3C.pdf"))
venn <- venn_analysis(peaks[ EBNA3C == 1],west, 0, "EBNA3A")
grid.draw(venn)
grid.newpage()
venn <- venn_analysis(peaks[ EBNA3C == 1 & Dnase == 1],west, 0, "EBNA3C\nand\nDHS")
grid.draw(venn)    
dev.off()


pdf(file = file.path(figs_dir,"figS1_all.pdf"))
venn <- venn_analysis(peaks[ EBNA3A == 1 | EBNA3B == 1 | EBNA3C == 1],west, 0, "EBNA3")
grid.draw(venn)
grid.newpage()
venn <- venn_analysis(peaks[ (EBNA3A == 1 | EBNA3B == 1 | EBNA3C == 1) & Dnase == 1],west, 0, "EBNA3\nand\nDHS")
grid.draw(venn)    
dev.off()
