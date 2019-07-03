load("clean.data/tree-abund.rda")
load("clean.data/shrub-abund.rda")
load("clean.data/herb-abund.rda")
load("clean.data/wood-phylo.rda")
load("clean.data/herb-phylo.rda")

hgt <- read.table("raw.data/hgt.txt")[[1]]

#1110-1180, 1190-1255, 1260-1338, 1344-1408, 1414-1496, 1504-1552, 1558-1616, 1622-1682, and 1687-1734
br <- c(1180, 1255, 1338, 1408, 1496, 1552, 1616, 1682)
#________________________________________

source("R/facet-tools.r")

compute.beta.c <- function(cdm, tree, q)
{
  res <- sapply(1:95, function(x) p.decomp(cdm[x:(x+1), ], tree, q))
  res[4, ]
}

d <- ta
tree <- wood.tree
q <- 2



system.time(c1 <- compute.beta.c.95(ha, herb.tree, 2))

system.time(c2 <- p.decomp.95(ha, herb.tree, 2)$c)





for (ii in 1:95) res <- p.decomp(ta[ii:(ii+1), ], wood.tree, 2)



library(picante)

null <- replicate(100, compute.beta.c(ta, tipShuffle(wood.tree), 0))
beta.ses <- data.frame(tc0 = (compute.beta.c(ta, wood.tree, 0) - rowMeans(null)) / apply(null, 1, sd))

null <- replicate(100, compute.beta.c(ta, tipShuffle(wood.tree), 1))
beta.ses$tc1 <- (compute.beta.c(ta, wood.tree, 1) - rowMeans(null)) / apply(null, 1, sd)

null <- replicate(100, compute.beta.c(ta, tipShuffle(wood.tree), 2))
beta.ses$tc2 <- (compute.beta.c(ta, wood.tree, 2) - rowMeans(null)) / apply(null, 1, sd)

null <- replicate(100, compute.beta.c(sa, tipShuffle(wood.tree), 0))
beta.ses$sc0 <- (compute.beta.c(sa, wood.tree, 0) - rowMeans(null)) / apply(null, 1, sd)

null <- replicate(100, compute.beta.c(sa, tipShuffle(wood.tree), 1))
beta.ses$sc1 <- (compute.beta.c(sa, wood.tree, 1) - rowMeans(null)) / apply(null, 1, sd)

null <- replicate(100, compute.beta.c(sa, tipShuffle(wood.tree), 2))
beta.ses$sc2 <- (compute.beta.c(sa, wood.tree, 2) - rowMeans(null)) / apply(null, 1, sd)

null <- replicate(100, compute.beta.c(ha, tipShuffle(herb.tree), 0))
beta.ses$hc0 <- (compute.beta.c(ha, herb.tree, 0) - rowMeans(null)) / apply(null, 1, sd)

null <- replicate(100, compute.beta.c(ha, tipShuffle(herb.tree), 1))
beta.ses$hc1 <- (compute.beta.c(ha, herb.tree, 1) - rowMeans(null)) / apply(null, 1, sd)

null <- replicate(100, compute.beta.c(ha, tipShuffle(herb.tree), 2))
beta.ses$hc2 <- (compute.beta.c(ha, herb.tree, 2) - rowMeans(null)) / apply(null, 1, sd)

save(beta.ses, file = "clean.data/pw-beta.rda")

source("R/plot.ses.r")

op <- par(mfcol = c(3, 3), mar = c(4, 4, 0.5, 0.5))
plot.ses(beta.ses$tc0, hgt[1:95], col = "forestgreen", lab = expression(PC[0]))
plot.ses(beta.ses$tc1, hgt[1:95], col = "forestgreen", lab = expression(PC[1]))
plot.ses(beta.ses$tc2, hgt[1:95], col = "forestgreen", lab = expression(PC[2]))

plot.ses(beta.ses$sc0, hgt[1:95], col = "skyblue", lab = expression(PC[0]))
plot.ses(beta.ses$sc1, hgt[1:95], col = "skyblue", lab = expression(PC[1]))
plot.ses(beta.ses$sc2, hgt[1:95], col = "skyblue", lab = expression(PC[2]))

plot.ses(beta.ses$hc0, hgt[1:95], col = "tomato", lab = expression(PC[0]))
plot.ses(beta.ses$hc1, hgt[1:95], col = "tomato", lab = expression(PC[1]))
plot.ses(beta.ses$hc2, hgt[1:95], col = "tomato", lab = expression(PC[2]))
par(op)


sapply(beta.ses, function(x) wilcox.test(x)$p.value)





beta <- data.frame(tc0 = compute.beta.c(ta, wood.tree, 0))
beta$tc1 <- compute.beta.c(ta, wood.tree, 1)
beta$tc2 <- compute.beta.c(ta, wood.tree, 2)
beta$sc0 <- compute.beta.c(sa, wood.tree, 0)
beta$sc1 <- compute.beta.c(sa, wood.tree, 1)
beta$sc2 <- compute.beta.c(sa, wood.tree, 2)
beta$hc0 <- compute.beta.c(ha, herb.tree, 0)
beta$hc1 <- compute.beta.c(ha, herb.tree, 1)
beta$hc2 <- compute.beta.c(ha, herb.tree, 2)

op <- par(mfcol = c(3, 3), mar = c(4, 4, 0.5, 0.5))
plot.ses(beta$tc0, hgt[1:95], col = "forestgreen", lab = expression(PC[0]))
plot.ses(beta$tc1, hgt[1:95], col = "forestgreen", lab = expression(PC[1]))
plot.ses(beta$tc2, hgt[1:95], col = "forestgreen", lab = expression(PC[2]))
abline(v = br, lty = 2, col = "grey")

plot.ses(beta$sc0, hgt[1:95], col = "skyblue", lab = expression(PC[0]))
plot.ses(beta$sc1, hgt[1:95], col = "skyblue", lab = expression(PC[1]))
plot.ses(beta$sc2, hgt[1:95], col = "skyblue", lab = expression(PC[2]))
abline(v = br, lty = 2, col = "grey")

plot.ses(beta$hc0, hgt[1:95], col = "tomato", lab = expression(PC[0]))
plot.ses(beta$hc1, hgt[1:95], col = "tomato", lab = expression(PC[1]))
plot.ses(beta$hc2, hgt[1:95], col = "tomato", lab = expression(PC[2]))
abline(v = br, lty = 2, col = "grey")
par(op)