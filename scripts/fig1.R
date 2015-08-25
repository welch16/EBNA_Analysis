
rm(list = ls())

library(VennDiagram)
library(data.table)

load(file = "data/RData/unified_lists_wProbs.RData") ## unified_lists

overlaps <- unified_lists$overlaps

area1 <- nrow(overlaps[EBNA2 == 1 & EBNA3A == 0 & EBNA3B == 0 & EBNA3C == 0])
area2 <- nrow(overlaps[EBNA2 == 0 & EBNA3A == 1 & EBNA3B == 0 & EBNA3C == 0])
area3 <- nrow(overlaps[EBNA2 == 0 & EBNA3A == 0 & EBNA3B == 1 & EBNA3C == 0])
area4 <- nrow(overlaps[EBNA2 == 0 & EBNA3A == 0 & EBNA3B == 0 & EBNA3C == 1])

n12 <- nrow(overlaps[EBNA2 == 1 & EBNA3A == 1 & EBNA3B == 0 & EBNA3C == 0])
n13 <- nrow(overlaps[EBNA2 == 1 & EBNA3A == 0 & EBNA3B == 1 & EBNA3C == 0])
n14 <- nrow(overlaps[EBNA2 == 1 & EBNA3A == 0 & EBNA3B == 0 & EBNA3C == 1])
n23 <- nrow(overlaps[EBNA2 == 0 & EBNA3A == 1 & EBNA3B == 1 & EBNA3C == 0])
n24 <- nrow(overlaps[EBNA2 == 0 & EBNA3A == 1 & EBNA3B == 0 & EBNA3C == 1])
n34 <- nrow(overlaps[EBNA2 == 0 & EBNA3A == 0 & EBNA3B == 1 & EBNA3C == 1])

n123 <- nrow(overlaps[EBNA2 == 1 & EBNA3A == 1 & EBNA3B == 1 & EBNA3C == 0])
n124 <- nrow(overlaps[EBNA2 == 1 & EBNA3A == 1 & EBNA3B == 0 & EBNA3C == 1])
n134 <- nrow(overlaps[EBNA2 == 1 & EBNA3A == 0 & EBNA3B == 1 & EBNA3C == 1])
n234 <- nrow(overlaps[EBNA2 == 0 & EBNA3A == 1 & EBNA3B == 1 & EBNA3C == 1])

n1234 <- nrow(overlaps[EBNA2 == 1 & EBNA3A == 1 & EBNA3B == 1 & EBNA3C == 1])

  ##   a6 <- n1234
  ##   a12 <- n123 - a6
  ##   a11 <- n124 - a6
  ##   a5 <- n134 - a6
  ##   a7 <- n234 - a6
  ##   a15 <- n12 - a6 - a11 - a12
  ##   a4 <- n13 - a6 - a5 - a12
  ##   a10 <- n14 - a6 - a5 - a11
  ##   a13 <- n23 - a6 - a7 - a12
  ##   a8 <- n24 - a6 - a7 - a11
  ##   a2 <- n34 - a6 - a5 - a7
  ##   a9 <- area1 - a4 - a5 - a6 - a10 - a11 - a12 - a15
  ##   a14 <- area2 - a6 - a7 - a8 - a11 - a12 - a13 - a15
  ##   a1 <- area3 - a2 - a4 - a5 - a6 - a7 - a12 - a13
  ##   a3 <- area4 - a2 - a5 - a6 - a7 - a8 - a10 - a11


  ## areas <- c(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, 
  ##       a12, a13, a14)

draw.quad.venn <- function (area1, area2, area3, area4, n12, n13, n14, n23, n24, 
    n34, n123, n124, n134, n234, n1234, category = rep("", 4), 
    lwd = rep(2, 4), lty = rep("solid", 4), col = rep("black", 
        4), fill = NULL, alpha = rep(0.5, 4), label.col = rep("black", 
        15), cex = rep(1, 15), fontface = rep("plain", 15), fontfamily = rep("serif", 
        15), cat.pos = c(-15, 15, 0, 0), cat.dist = c(0.22, 0.22, 
        0.11, 0.11), cat.col = rep("black", 4), cat.cex = rep(1, 
        4), cat.fontface = rep("plain", 4), cat.fontfamily = rep("serif", 
        4), cat.just = rep(list(c(0.5, 0.5)), 4), rotation.degree = 0, 
    rotation.centre = c(0.5, 0.5), ind = TRUE, ...) 
{
    if (length(category) == 1) {
        cat <- rep(category, 4)
    }
    else if (length(category) != 4) {
        stop("Unexpected parameter length for 'category'")
    }
    if (length(lwd) == 1) {
        lwd <- rep(lwd, 4)
    }
    else if (length(lwd) != 4) {
        stop("Unexpected parameter length for 'lwd'")
    }
    if (length(lty) == 1) {
        lty <- rep(lty, 4)
    }
    else if (length(lty) != 4) {
        stop("Unexpected parameter length for 'lty'")
    }
    if (length(col) == 1) {
        col <- rep(col, 4)
    }
    else if (length(col) != 4) {
        stop("Unexpected parameter length for 'col'")
    }
    if (length(label.col) == 1) {
        label.col <- rep(label.col, 15)
    }
    else if (length(label.col) != 15) {
        stop("Unexpected parameter length for 'label.col'")
    }
    if (length(cex) == 1) {
        cex <- rep(cex, 15)
    }
    else if (length(cex) != 15) {
        stop("Unexpected parameter length for 'cex'")
    }
    if (length(fontface) == 1) {
        fontface <- rep(fontface, 15)
    }
    else if (length(fontface) != 15) {
        stop("Unexpected parameter length for 'fontface'")
    }
    if (length(fontfamily) == 1) {
        fontfamily <- rep(fontfamily, 15)
    }
    else if (length(fontfamily) != 15) {
        stop("Unexpected parameter length for 'fontfamily'")
    }
    if (length(fill) == 1) {
        fill <- rep(fill, 4)
    }
    else if (length(fill) != 4 & length(fill) != 0) {
        stop("Unexpected parameter length for 'fill'")
    }
    if (length(alpha) == 1) {
        alpha <- rep(alpha, 4)
    }
    else if (length(alpha) != 4 & length(alpha) != 0) {
        stop("Unexpected parameter length for 'alpha'")
    }
    if (length(cat.pos) == 1) {
        cat.pos <- rep(cat.pos, 4)
    }
    else if (length(cat.pos) != 4) {
        stop("Unexpected parameter length for 'cat.pos'")
    }
    if (length(cat.dist) == 1) {
        cat.dist <- rep(cat.dist, 4)
    }
    else if (length(cat.dist) != 4) {
        stop("Unexpected parameter length for 'cat.dist'")
    }
    if (length(cat.col) == 1) {
        cat.col <- rep(cat.col, 4)
    }
    else if (length(cat.col) != 4) {
        stop("Unexpected parameter length for 'cat.col'")
    }
    if (length(cat.cex) == 1) {
        cat.cex <- rep(cat.cex, 4)
    }
    else if (length(cat.cex) != 4) {
        stop("Unexpected parameter length for 'cat.cex'")
    }
    if (length(cat.fontface) == 1) {
        cat.fontface <- rep(cat.fontface, 4)
    }
    else if (length(cat.fontface) != 4) {
        stop("Unexpected parameter length for 'cat.fontface'")
    }
    if (length(cat.fontfamily) == 1) {
        cat.fontfamily <- rep(cat.fontfamily, 4)
    }
    else if (length(cat.fontfamily) != 4) {
        stop("Unexpected parameter length for 'cat.fontfamily'")
    }
    if (!(class(cat.just) == "list" & length(cat.just) == 4 & 
        length(cat.just[[1]]) == 2 & length(cat.just[[2]]) == 
        2 & length(cat.just[[3]]) == 2 & length(cat.just[[4]]) == 
        2)) {
        stop("Unexpected parameter format for 'cat.just'")
    }
    cat.pos <- cat.pos + rotation.degree

    grob.list <- gList()
    ellipse.positions <- matrix(nrow = 4, ncol = 7)
    colnames(ellipse.positions) <- c("x", "y", "a", "b", "rotation", 
        "fill.mapping", "line.mapping")
    ellipse.positions[1, ] <- c(0.65, 0.47, 0.35, 0.2, 45, 2, 
        4)
    ellipse.positions[2, ] <- c(0.35, 0.47, 0.35, 0.2, 135, 1, 
        1)
    ellipse.positions[3, ] <- c(0.5, 0.57, 0.33, 0.15, 45, 4, 
        3)
    ellipse.positions[4, ] <- c(0.5, 0.57, 0.35, 0.15, 135, 3, 
        2)
    for (i in 1:4) {
        grob.list <- gList(grob.list, VennDiagram::ellipse(x = ellipse.positions[i, 
            "x"], y = ellipse.positions[i, "y"], a = ellipse.positions[i, 
            "a"], b = ellipse.positions[i, "b"], rotation = ellipse.positions[i, 
            "rotation"], gp = gpar(lty = 0, fill = fill[ellipse.positions[i, 
            "fill.mapping"]], alpha = alpha[ellipse.positions[i, 
            "fill.mapping"]])))
    }
    for (i in 1:4) {
        grob.list <- gList(grob.list, ellipse(x = ellipse.positions[i, 
            "x"], y = ellipse.positions[i, "y"], a = ellipse.positions[i, 
            "a"], b = ellipse.positions[i, "b"], rotation = ellipse.positions[i, 
            "rotation"], gp = gpar(lwd = lwd[ellipse.positions[i, 
            "line.mapping"]], lty = lty[ellipse.positions[i, 
            "line.mapping"]], col = col[ellipse.positions[i, 
            "line.mapping"]], fill = "transparent")))
    }
    label.matrix <- matrix(nrow = 15, ncol = 3)
    colnames(label.matrix) <- c("label", "x", "y")

    label.matrix[1, ] <- c(area3, 0.35, 0.77)
    label.matrix[2, ] <- c(n34, 0.5, 0.69)
    label.matrix[3, ] <- c(area4, 0.65, 0.77)
    label.matrix[4, ] <- c(n13, 0.31, 0.67)
    label.matrix[5, ] <- c(n134, 0.4, 0.58)
    label.matrix[6, ] <- c(n1234, 0.5, 0.47)
    label.matrix[7, ] <- c(n234, 0.6, 0.58)
    label.matrix[8, ] <- c(n24, 0.69, 0.67)
    label.matrix[9, ] <- c(area1, 0.18, 0.58)
    label.matrix[10, ] <- c(n14, 0.32, 0.42)
    label.matrix[11, ] <- c(n124, 0.425, 0.38)
    label.matrix[12, ] <- c(n234, 0.575, 0.38)
    label.matrix[13, ] <- c(n23, 0.68, 0.42)
    label.matrix[14, ] <- c(area2, 0.82, 0.58)
    label.matrix[15, ] <- c(n12, 0.5, 0.28)
    
    for (i in 1:nrow(label.matrix)) {
        grob.list <- gList(grob.list, textGrob(label = label.matrix[i, 
            "label"], x = label.matrix[i, "x"], y = label.matrix[i, 
            "y"], gp = gpar(col = label.col[i], cex = cex[i], 
            fontface = fontface[i], fontfamily = fontfamily[i])))
    }
    cat.pos.x <- c(0.18, 0.82, 0.35, 0.65)
    cat.pos.y <- c(0.58, 0.58, 0.77, 0.77)
    for (i in 1:4) {
        this.cat.pos <- find.cat.pos(x = cat.pos.x[i], y = cat.pos.y[i], 
            pos = cat.pos[i], dist = cat.dist[i])
        grob.list <- gList(grob.list, textGrob(label = category[i], 
            x = this.cat.pos$x, y = this.cat.pos$y, just = cat.just[[i]], 
            gp = gpar(col = cat.col[i], cex = cat.cex[i], fontface = cat.fontface[i], 
                fontfamily = cat.fontfamily[i])))
    }
    grob.list <- VennDiagram::adjust.venn(VennDiagram::rotate.venn.degrees(grob.list, 
        rotation.degree, rotation.centre[1], rotation.centre[2]), 
        ...)
    if (ind) {
        grid.draw(grob.list)
    }
    return(grob.list)
}

category <- c("EBNA2","EBNA3A","EBNA3B","EBNA3C")

venn <- draw.quad.venn(area1 = area1,
                       area2 = area2,
                       area3 = area3,
                       area4 = area4,
                       n12 = n12,
                       n13 = n13,
                       n14 = n14,
                       n23 = n23,
                       n24 = n24,
                       n34 = n34,
                       n123 = n123,
                       n124 = n124,
                       n134 = n134,
                       n234 = n234,
                       n1234 = n1234,
                       category = category,
                       cex = 1.8,
                       lwd = 3,
                       cat.cex = 1.8,
                       col = c("purple","blue","red","darkgreen"))
pdf(file = "figures/for_paper/fig1.pdf")
grid.draw(venn)
dev.off()


