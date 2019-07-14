library(picante)
library(GUniFrac)

load("clean.data/tree-abund.rda")
load("clean.data/shrub-abund.rda")
load("clean.data/herb-abund.rda")
load("clean.data/wood-phylo.rda")
load("clean.data/herb-phylo.rda")

hgt <- read.table("raw.data/hgt.txt")[[1]]

shrub.tree <- drop.tip(wood.tree, wood.tree$tip.label[!(wood.tree$tip.label %in% colnames(sa))])
wood.tree <- drop.tip(wood.tree, wood.tree$tip.label[!(wood.tree$tip.label %in% colnames(ta))])

#1110-1180, 1190-1255, 1260-1338, 1344-1408, 1414-1496, 1504-1552, 1558-1616, 1622-1682, and 1687-1734
br <- c(1180, 1255, 1338, 1408, 1496, 1552, 1616, 1682)

#________________________________________

uf.rand <- function(cdm, tree, ab)
{
  ntree <- tipShuffle(tree)
  if (ab){
    sapply(1:95, function(x) (GUniFrac(cdm[x:(x+1), ], ntree, alpha = 1)$unifracs)[, , "d_1"])[3, ]
  } else {
    sapply(1:95, function(x) (GUniFrac(cdm[x:(x+1), ], ntree, alpha = 1)$unifracs)[, , "d_UW"])[3, ]
  }
}

null <- replicate(100, uf.rand(ha, herb.tree, F))
beta.ses <- data.frame(hd0 = ((sapply(1:95, function(x) (GUniFrac(ha[x:(x+1), ],
                              herb.tree, alpha = 1)$unifracs)[, , "d_UW"])[3, ] - rowMeans(null)) / apply(null, 1, sd)))

null <- replicate(100, uf.rand(ha, herb.tree, T))
beta.ses$hd1 <- (sapply(1:95, function(x) (GUniFrac(ha[x:(x+1), ],
                 herb.tree, alpha = 1)$unifracs)[, , "d_1"])[3, ] - rowMeans(null)) / apply(null, 1, sd)



null <- replicate(100, uf.rand(ta, wood.tree, F))
beta.ses$td0 <- (sapply(1:95, function(x) (GUniFrac(ta[x:(x+1), ],
                 wood.tree, alpha = 1)$unifracs)[, , "d_UW"])[3, ] - rowMeans(null)) / apply(null, 1, sd)

null <- replicate(100, uf.rand(ta, wood.tree, T))
beta.ses$td1 <- (sapply(1:95, function(x) (GUniFrac(ta[x:(x+1), ],
                 wood.tree, alpha = 1)$unifracs)[, , "d_1"])[3, ] - rowMeans(null)) / apply(null, 1, sd)



null <- replicate(100, uf.rand(sa, shrub.tree, F))
beta.ses$sd0 <- (sapply(1:95, function(x) (GUniFrac(sa[x:(x+1), ],
                 shrub.tree, alpha = 1)$unifracs)[, , "d_UW"])[3, ] - rowMeans(null)) / apply(null, 1, sd)

null <- replicate(100, uf.rand(sa, shrub.tree, T))
beta.ses$sd1 <- (sapply(1:95, function(x) (GUniFrac(sa[x:(x+1), ],
                 shrub.tree, alpha = 1)$unifracs)[, , "d_1"])[3, ] - rowMeans(null)) / apply(null, 1, sd)

save(beta.ses, file = "clean.data/pw-beta-uf-tl.rda")




load("clean.data/pw-beta-uf-tl.rda")

source("R/plot.ses.r")

op <- par(mfcol = c(2, 3), mar = c(4, 4, 0.5, 0.5))
plot.ses(beta.ses$td0[c(T, F)], hgt[1:95][c(T, F)], col = "forestgreen", lab = expression(PC[0]))
plot.ses(beta.ses$td1[c(T, F)], hgt[1:95][c(T, F)], col = "forestgreen", lab = expression(PC[1]))

plot.ses(beta.ses$sd0[c(T, F)], hgt[1:95][c(T, F)], col = "skyblue", lab = expression(PC[0]))
plot.ses(beta.ses$sd1[c(T, F)], hgt[1:95][c(T, F)], col = "skyblue", lab = expression(PC[1]))

plot.ses(beta.ses$hd0[c(T, F)], hgt[1:95][c(T, F)], col = "tomato", lab = expression(PC[0]))
plot.ses(beta.ses$hd1[c(T, F)], hgt[1:95][c(T, F)], col = "tomato", lab = expression(PC[1]))
par(op)
