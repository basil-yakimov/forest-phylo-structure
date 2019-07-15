load("clean.data/tree-abund.rda")
load("clean.data/shrub-abund.rda")
load("clean.data/herb-abund.rda")
load("clean.data/wood-phylo.rda")
load("clean.data/herb-phylo.rda")

hgt <- read.table("raw.data/hgt.txt")[[1]]

#1110-1180, 1190-1255, 1260-1338, 1344-1408, 1414-1496, 1504-1552, 1558-1616, 1622-1682, and 1687-1734
br <- c(1180, 1255, 1338, 1408, 1496, 1552, 1616, 1682)
#________________________________________

#dpw <- sapply(1:95, function(x) comdist(sa[x:(x+1), ], wood.mat, T))

load("clean.data/pw-beta-tl.rda")

library(picante)
wood.mat <- cophenetic(wood.tree)
herb.mat <- cophenetic(herb.tree)

dnn.rand <- function(cdm, mat, ab)
{
  mat <- mat[row.names(mat) %in% colnames(cdm), row.names(mat) %in% colnames(cdm)]
  ord <- sample(nrow(mat))
  row.names(mat) <- row.names(mat)[ord]
  colnames(mat) <- colnames(mat)[ord]
  sapply(1:95, function(x) comdistnt(cdm[x:(x+1), ], mat[ord, ord], ab))
}

null <- replicate(1000, dnn.rand(ta, wood.mat, F))
beta.ses$tdnn0 <- (sapply(1:95, function(x) comdistnt(ta[x:(x+1), ], wood.mat, F)) - rowMeans(null)) / apply(null, 1, sd)

null <- replicate(1000, dnn.rand(ta, wood.mat, T))
beta.ses$tdnn1 <- (sapply(1:95, function(x) comdistnt(ta[x:(x+1), ], wood.mat, T)) - rowMeans(null)) / apply(null, 1, sd)

null <- replicate(1000, dnn.rand(sa, wood.mat, F))
beta.ses$sdnn0 <- (sapply(1:95, function(x) comdistnt(sa[x:(x+1), ], wood.mat, F)) - rowMeans(null)) / apply(null, 1, sd)

null <- replicate(1000, dnn.rand(sa, wood.mat, T))
beta.ses$sdnn1 <- (sapply(1:95, function(x) comdistnt(sa[x:(x+1), ], wood.mat, T)) - rowMeans(null)) / apply(null, 1, sd)

null <- replicate(1000, dnn.rand(ha, herb.mat, F))
beta.ses$hdnn0 <- (sapply(1:95, function(x) comdistnt(ha[x:(x+1), ], herb.mat, F)) - rowMeans(null)) / apply(null, 1, sd)

null <- replicate(1000, dnn.rand(ha, herb.mat, T))
beta.ses$hdnn1 <- (sapply(1:95, function(x) comdistnt(ha[x:(x+1), ], herb.mat, T)) - rowMeans(null)) / apply(null, 1, sd)


save(beta.ses, file = "clean.data/pw-beta-tl.rda")

source("R/plot.ses.r")

op <- par(mfcol = c(2, 3), mar = c(4, 4, 0.5, 0.5))
plot.ses(beta.ses$tdnn0[c(T, F)], hgt[1:95][c(T, F)], col = "forestgreen", lab = expression(PC[0]))
plot.ses(beta.ses$tdnn1[c(T, F)], hgt[1:95][c(T, F)], col = "forestgreen", lab = expression(PC[1]))

plot.ses(beta.ses$sdnn0[c(T, F)], hgt[1:95][c(T, F)], col = "skyblue", lab = expression(PC[0]))
plot.ses(beta.ses$sdnn1[c(T, F)], hgt[1:95][c(T, F)], col = "skyblue", lab = expression(PC[1]))

plot.ses(beta.ses$hdnn0[c(T, F)], hgt[1:95][c(T, F)], col = "tomato", lab = expression(PC[0]))
plot.ses(beta.ses$hdnn1[c(T, F)], hgt[1:95][c(T, F)], col = "tomato", lab = expression(PC[1]))
par(op)


sapply(beta.ses, function(x) wilcox.test(x)$p.value)