library(abind)
library(picante)

load("clean.data/tree-abund.rda")
load("clean.data/shrub-abund.rda")
load("clean.data/herb-abund.rda")
load("clean.data/wood-phylo.rda")
load("clean.data/herb-phylo.rda")

hgt <- read.table("raw.data/hgt.txt")[[1]]



beta <- data.frame(tp2 = as.vector(phylosor(ta, wood.tree)))
null <- phylosor.rnd(ta, wood.tree, cstSor=TRUE, null.model = "taxa.labels", runs = 100)
beta$tp2.null <- rowMeans(matrix(unlist(null), ncol = 100, byrow = TRUE))
beta$tp2.ses <- (beta$tp2 - beta$tp2.null)/apply(matrix(unlist(null), ncol = 100, byrow = TRUE), 1, sd)

beta$sp2 <- as.vector(phylosor(sa, wood.tree))
null <- phylosor.rnd(sa, wood.tree, cstSor=TRUE, null.model = "taxa.labels", runs = 100)
beta$sp2.null <- rowMeans(matrix(unlist(null), ncol = 100, byrow = TRUE))
beta$sp2.ses <- (beta$sp2 - beta$sp2.null)/apply(matrix(unlist(null), ncol = 100, byrow = TRUE), 1, sd)

beta$hp2 <- as.vector(phylosor(ha, herb.tree))
null <- phylosor.rnd(ha, herb.tree, cstSor=TRUE, null.model = "taxa.labels", runs = 100)
beta$hp2.null <- rowMeans(matrix(unlist(null), ncol = 100, byrow = TRUE))
beta$hp2.ses <- (beta$hp2 - beta$hp2.null)/apply(matrix(unlist(null), ncol = 100, byrow = TRUE), 1, sd)

beta$dh <- dist(hgt)
beta$mh <- as.dist(sapply(hgt, function(x) sapply(hgt, function(z) (x+z)/2)))

save(beta, file = "clean.data/beta-using-picante.rda")


#_________________________________________________________________


load("clean.data/beta-using-picante.rda")

source("R/plot.ses.r")

png("figures/beta-phylo-picante.png", width = 3000, height = 2000)
op <- par(mfcol = c(2, 3), mar = c(4, 4, 0.5, 0.5), cex = 2)

plot.ses(beta$tp2, beta$dh, col = rgb(0.25, 0.75, 0.25, 0.3), lab = expression("phylo - " * beta[2]), xlab = expression(Delta * "hgt"))
plot.ses(beta$tp2, beta$mh, col = rgb(0.25, 0.75, 0.25, 0.3), lab = expression("phylo - " * beta[2]), xlab = "mean hgt")

plot.ses(beta$sp2, beta$dh, col = rgb(0.25, 0.25, 0.75, 0.3), lab = expression("phylo - " * beta[2]), xlab = expression(Delta * "hgt"))
plot.ses(beta$sp2, beta$mh, col = rgb(0.25, 0.25, 0.75, 0.3), lab = expression("phylo - " * beta[2]), xlab = "mean hgt")

plot.ses(beta$hp2, beta$dh, col = rgb(0.75, 0.25, 0.25, 0.3), lab = expression("phylo - " * beta[2]), xlab = expression(Delta * "hgt"))
plot.ses(beta$hp2, beta$mh, col = rgb(0.75, 0.25, 0.25, 0.3), lab = expression("phylo - " * beta[2]), xlab = "mean hgt")

par(op)
dev.off()


