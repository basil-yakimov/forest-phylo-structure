load("clean.data/tree-abund.rda")
load("clean.data/shrub-abund.rda")
load("clean.data/herb-abund.rda")
load("clean.data/wood-phylo.rda")
load("clean.data/herb-phylo.rda")

hgt <- read.table("raw.data/hgt.txt")[[1]]

#1110-1180,1190-1255, 1260-1338, 1344-1408, 1414-1496, 1504-1552, 1558-1616, 1622-1682, and 1687-1734
subs <- cut(hgt, breaks = c(0, 1185, 1257, 1340, 1410, 1500, 1554, 1620, 1685, 3000),
            labels = paste0("sub", 1:9))

table(subs)

#________________________________________

source("R/facet-tools.r")

compute.beta.c <- function(cdm, tree, q)
{
  res <- tapply(1:96, subs, function(x) p.decomp(cdm[x, ], tree, q))
  simplify2array(res)[4, ]
}

compute.beta.c(ha, herb.tree, 2)

library(picante)

null <- replicate(1000, compute.beta.c(randomizeMatrix(ta, null.model = "independentswap"), wood.tree, 0))
beta.ses <- data.frame(tc0 = (compute.beta.c(ta, wood.tree, 0) - rowMeans(null)) / apply(null, 1, sd))

null <- replicate(1000, compute.beta.c(randomizeMatrix(ta, null.model = "independentswap"), wood.tree, 1))
beta.ses$tc1 <- (compute.beta.c(ta, wood.tree, 1) - rowMeans(null)) / apply(null, 1, sd)

null <- replicate(1000, compute.beta.c(randomizeMatrix(ta, null.model = "independentswap"), wood.tree, 2))
beta.ses$tc2 <- (compute.beta.c(ta, wood.tree, 2) - rowMeans(null)) / apply(null, 1, sd)

null <- replicate(1000, compute.beta.c(randomizeMatrix(sa, null.model = "independentswap"), wood.tree, 0))
beta.ses$sc0 <- (compute.beta.c(sa, wood.tree, 0) - rowMeans(null)) / apply(null, 1, sd)

null <- replicate(1000, compute.beta.c(randomizeMatrix(sa, null.model = "independentswap"), wood.tree, 1))
beta.ses$sc1 <- (compute.beta.c(sa, wood.tree, 1) - rowMeans(null)) / apply(null, 1, sd)

null <- replicate(1000, compute.beta.c(randomizeMatrix(sa, null.model = "independentswap"), wood.tree, 2))
beta.ses$sc2 <- (compute.beta.c(sa, wood.tree, 2) - rowMeans(null)) / apply(null, 1, sd)

null <- replicate(1000, compute.beta.c(randomizeMatrix(ha, null.model = "independentswap"), herb.tree, 0))
beta.ses$hc0 <- (compute.beta.c(ha, herb.tree, 0) - rowMeans(null)) / apply(null, 1, sd)

null <- replicate(1000, compute.beta.c(randomizeMatrix(ha, null.model = "independentswap"), herb.tree, 1))
beta.ses$hc1 <- (compute.beta.c(ha, herb.tree, 1) - rowMeans(null)) / apply(null, 1, sd)

null <- replicate(1000, compute.beta.c(randomizeMatrix(ha, null.model = "independentswap"), herb.tree, 2))
beta.ses$hc2 <- (compute.beta.c(ha, herb.tree, 2) - rowMeans(null)) / apply(null, 1, sd)



source("R/plot.ses.r")
mhgt <- tapply(hgt, subs, mean)

op <- par(mfcol = c(3, 3), mar = c(4, 4, 0.5, 0.5))
plot.ses(beta.ses$tc0, mhgt, col = "forestgreen", lab = expression(PC[0]))
plot.ses(beta.ses$tc1, mhgt, col = "forestgreen", lab = expression(PC[1]))
plot.ses(beta.ses$tc2, mhgt, col = "forestgreen", lab = expression(PC[2]))

plot.ses(beta.ses$sc0, mhgt, col = "skyblue", lab = expression(PC[0]))
plot.ses(beta.ses$sc1, mhgt, col = "skyblue", lab = expression(PC[1]))
plot.ses(beta.ses$sc2, mhgt, col = "skyblue", lab = expression(PC[2]))

plot.ses(beta.ses$hc0, mhgt, col = "tomato", lab = expression(PC[0]))
plot.ses(beta.ses$hc1, mhgt, col = "tomato", lab = expression(PC[1]))
plot.ses(beta.ses$hc2, mhgt, col = "tomato", lab = expression(PC[2]))
par(op)


sapply(beta.ses, function(x) wilcox.test(x)$p.value)
