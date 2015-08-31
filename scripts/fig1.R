
rm(list = ls())

library(VennDiagram)
library(data.table)

load(file = "data/RData/unified_lists_wProbs.RData") ## unified_lists
load("data/RData/factor_overlaps.RData")

source("R/venn_functions.R")

overlaps <- unified_lists$overlaps

category <- c("EBNA2","EBNA3A","EBNA3B","EBNA3C")

venn <- draw.quad.venn(area1 = nrow(overlaps[EBNA2 == 1 & EBNA3A == 0 & EBNA3B == 0 & EBNA3C == 0]),
                       area2 = nrow(overlaps[EBNA2 == 0 & EBNA3A == 1 & EBNA3B == 0 & EBNA3C == 0]),
                       area3 = nrow(overlaps[EBNA2 == 0 & EBNA3A == 0 & EBNA3B == 1 & EBNA3C == 0]),
                       area4 = nrow(overlaps[EBNA2 == 0 & EBNA3A == 0 & EBNA3B == 0 & EBNA3C == 1]),
                       n12 = nrow(overlaps[EBNA2 == 1 & EBNA3A == 1 & EBNA3B == 0 & EBNA3C == 0]),
                       n13 = nrow(overlaps[EBNA2 == 1 & EBNA3A == 0 & EBNA3B == 1 & EBNA3C == 0]),
                       n14 = nrow(overlaps[EBNA2 == 1 & EBNA3A == 0 & EBNA3B == 0 & EBNA3C == 1]),
                       n23 = nrow(overlaps[EBNA2 == 0 & EBNA3A == 1 & EBNA3B == 1 & EBNA3C == 0]),
                       n24 = nrow(overlaps[EBNA2 == 0 & EBNA3A == 1 & EBNA3B == 0 & EBNA3C == 1]),
                       n34 = nrow(overlaps[EBNA2 == 0 & EBNA3A == 0 & EBNA3B == 1 & EBNA3C == 1]),
                       n123 = nrow(overlaps[EBNA2 == 1 & EBNA3A == 1 & EBNA3B == 1 & EBNA3C == 0]),
                       n124 = nrow(overlaps[EBNA2 == 1 & EBNA3A == 1 & EBNA3B == 0 & EBNA3C == 1]),
                       n134 = nrow(overlaps[EBNA2 == 1 & EBNA3A == 0 & EBNA3B == 1 & EBNA3C == 1]),
                       n234 = nrow(overlaps[EBNA2 == 0 & EBNA3A == 1 & EBNA3B == 1 & EBNA3C == 1]),
                       n1234 = nrow(overlaps[EBNA2 == 1 & EBNA3A == 1 & EBNA3B == 1 & EBNA3C == 1]),
                       category = category,
                       cex = 1.8,
                       lwd = 3,
                       cat.cex = 1.8,
                       col = c("purple","blue","red","darkgreen"))

overlaps <- cbind(overlaps,factor_overlaps)[Dnase == 1]

venn_dnase <- draw.quad.venn(area1 = nrow(overlaps[EBNA2 == 1 & EBNA3A == 0 & EBNA3B == 0 & EBNA3C == 0]),
                       area2 = nrow(overlaps[EBNA2 == 0 & EBNA3A == 1 & EBNA3B == 0 & EBNA3C == 0]),
                       area3 = nrow(overlaps[EBNA2 == 0 & EBNA3A == 0 & EBNA3B == 1 & EBNA3C == 0]),
                       area4 = nrow(overlaps[EBNA2 == 0 & EBNA3A == 0 & EBNA3B == 0 & EBNA3C == 1]),
                       n12 = nrow(overlaps[EBNA2 == 1 & EBNA3A == 1 & EBNA3B == 0 & EBNA3C == 0]),
                       n13 = nrow(overlaps[EBNA2 == 1 & EBNA3A == 0 & EBNA3B == 1 & EBNA3C == 0]),
                       n14 = nrow(overlaps[EBNA2 == 1 & EBNA3A == 0 & EBNA3B == 0 & EBNA3C == 1]),
                       n23 = nrow(overlaps[EBNA2 == 0 & EBNA3A == 1 & EBNA3B == 1 & EBNA3C == 0]),
                       n24 = nrow(overlaps[EBNA2 == 0 & EBNA3A == 1 & EBNA3B == 0 & EBNA3C == 1]),
                       n34 = nrow(overlaps[EBNA2 == 0 & EBNA3A == 0 & EBNA3B == 1 & EBNA3C == 1]),
                       n123 = nrow(overlaps[EBNA2 == 1 & EBNA3A == 1 & EBNA3B == 1 & EBNA3C == 0]),
                       n124 = nrow(overlaps[EBNA2 == 1 & EBNA3A == 1 & EBNA3B == 0 & EBNA3C == 1]),
                       n134 = nrow(overlaps[EBNA2 == 1 & EBNA3A == 0 & EBNA3B == 1 & EBNA3C == 1]),
                       n234 = nrow(overlaps[EBNA2 == 0 & EBNA3A == 1 & EBNA3B == 1 & EBNA3C == 1]),
                       n1234 = nrow(overlaps[EBNA2 == 1 & EBNA3A == 1 & EBNA3B == 1 & EBNA3C == 1]),
                       category = category,
                       cex = 1.8,
                       lwd = 3,
                       cat.cex = 1.8,
                       col = c("purple","blue","red","darkgreen"))

pdf(file = "figures/for_paper/fig1.pdf")
grid.draw(venn)
grid.newpage()
grid.draw(venn_dnase)
dev.off()




