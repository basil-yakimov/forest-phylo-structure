load("clean.data/tree-abund.rda")
load("clean.data/shrub-abund.rda")
load("clean.data/herb-abund.rda")
load("clean.data/wood-phylo.rda")
load("clean.data/herb-phylo.rda")

hgt <- read.table("raw.data/hgt.txt")[[1]]

#1110-1180, 1190-1255, 1260-1338, 1344-1408, 1414-1496, 1504-1552, 1558-1616, 1622-1682, and 1687-1734
br <- c(1180, 1255, 1338, 1408, 1496, 1552, 1616, 1682)
#________________________________________

library(picante)
source("R/facet-tools.r")

tp0 <- p.dist(ta, wood.tree, 0)
tp1 <- p.dist(ta, wood.tree, 1)
tp2 <- p.dist(ta, wood.tree, 2)

beta <- data.frame(tp0 = as.vector(p.dist(ta, wood.tree, 0)))
null <- replicate(100, p.dist(ta, tipShuffle(wood.tree), 0))
beta$tp0.null <- rowMeans(null)
beta$tp0.ses <- (beta$tp0 - beta$tp0.null)/apply(null, 1, sd)

beta$tp1 = as.vector(p.dist(ta, wood.tree, 1))
null <- replicate(100, p.dist(ta, tipShuffle(wood.tree), 1))
beta$tp1.null <- rowMeans(null)
beta$tp1.ses <- (beta$tp1 - beta$tp1.null)/apply(null, 1, sd)

beta$tp2 = as.vector(p.dist(ta, wood.tree, 2))
null <- replicate(100, p.dist(ta, tipShuffle(wood.tree), 2))
beta$tp2.null <- rowMeans(null)
beta$tp2.ses <- (beta$tp2 - beta$tp2.null)/apply(null, 1, sd)

beta$sp0 = as.vector(p.dist(sa, wood.tree, 0))
null <- replicate(100, p.dist(sa, tipShuffle(wood.tree), 0))
beta$sp0.null <- rowMeans(null)
beta$sp0.ses <- (beta$sp0 - beta$sp0.null)/apply(null, 1, sd)

beta$sp1 = as.vector(p.dist(sa, wood.tree, 1))
null <- replicate(100, p.dist(sa, tipShuffle(wood.tree), 1))
beta$sp1.null <- rowMeans(null)
beta$sp1.ses <- (beta$sp1 - beta$sp1.null)/apply(null, 1, sd)

beta$sp2 = as.vector(p.dist(sa, wood.tree, 2))
null <- replicate(100, p.dist(sa, tipShuffle(wood.tree), 2))
beta$sp2.null <- rowMeans(null)
beta$sp2.ses <- (beta$sp2 - beta$sp2.null)/apply(null, 1, sd)


beta$hp0 = as.vector(p.dist(ha, herb.tree, 0))
null <- replicate(100, p.dist(ha, tipShuffle(herb.tree), 0))
beta$hp0.null <- rowMeans(null)
beta$hp0.ses <- (beta$hp0 - beta$hp0.null)/apply(null, 1, sd)

beta$hp1 = as.vector(p.dist(ha, herb.tree, 1))
null <- replicate(100, p.dist(ha, tipShuffle(herb.tree), 1))
beta$hp1.null <- rowMeans(null)
beta$hp1.ses <- (beta$hp1 - beta$hp1.null)/apply(null, 1, sd)

beta$hp2 = as.vector(p.dist(ha, herb.tree, 2))
null <- replicate(100, p.dist(ha, tipShuffle(herb.tree), 2))
beta$hp2.null <- rowMeans(null)
beta$hp2.ses <- (beta$hp2 - beta$hp2.null)/apply(null, 1, sd)


beta$tt0 <- t.dist(ta, 0)
beta$tt1 <- t.dist(ta, 1)
beta$tt2 <- t.dist(ta, 2)
beta$st0 <- t.dist(sa, 0)
beta$st1 <- t.dist(sa, 1)
beta$st2 <- t.dist(sa, 2)
beta$ht0 <- t.dist(ha, 0)
beta$ht1 <- t.dist(ha, 1)
beta$ht2 <- t.dist(ha, 2)

beta$dh <- dist(hgt)
beta$mh <- as.dist(sapply(hgt, function(x) sapply(hgt, function(z) (x+z)/2)))

save(beta, file = "clean.data/beta.full.rda")