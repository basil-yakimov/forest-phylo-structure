load("clean.data/tree-abund.rda")
load("clean.data/shrub-abund.rda")
load("clean.data/herb-abund.rda")
load("clean.data/wood-phylo.rda")
load("clean.data/herb-phylo.rda")

hgt <- read.table("raw.data/hgt.txt")[[1]]

#________________________________________

library(picante)
wood.mat <- cophenetic(wood.tree)
herb.mat <- cophenetic(herb.tree)

dpw.rand <- function(cdm, mat, ab)
{
  mat <- mat[row.names(mat) %in% colnames(cdm), row.names(mat) %in% colnames(cdm)]
  ord <- sample(nrow(mat))
  row.names(mat) <- row.names(mat)[ord]
  colnames(mat) <- colnames(mat)[ord]
  as.vector(comdist(cdm, mat[ord, ord], ab))
}



beta.dpw <- data.frame(td0 = as.vector(comdist(ta, wood.mat, F)))
null <- replicate(100, dpw.rand(ta, wood.mat, F))
beta.dpw$td0.null <- rowMeans(null)
beta.dpw$td0.ses <- (beta.dpw$td0 - beta.dpw$td0.null)/apply(null, 1, sd)

beta.dpw$td1 <- as.vector(comdist(ta, wood.mat, T))
null <- replicate(100, dpw.rand(ta, wood.mat, T))
beta.dpw$td1.null <- rowMeans(null)
beta.dpw$td1.ses <- (beta.dpw$td1 - beta.dpw$td1.null)/apply(null, 1, sd)


beta.dpw$sd0 <- as.vector(comdist(sa, wood.mat, F))
null <- replicate(100, dpw.rand(sa, wood.mat, F))
beta.dpw$sd0.null <- rowMeans(null)
beta.dpw$sd0.ses <- (beta.dpw$sd0 - beta.dpw$sd0.null)/apply(null, 1, sd)

beta.dpw$sd1 <- as.vector(comdist(sa, wood.mat, T))
null <- replicate(100, dpw.rand(sa, wood.mat, T))
beta.dpw$sd1.null <- rowMeans(null)
beta.dpw$sd1.ses <- (beta.dpw$sd1 - beta.dpw$sd1.null)/apply(null, 1, sd)


beta.dpw$hd0 <- as.vector(comdist(ha, herb.mat, F))
null <- replicate(100, dpw.rand(ha, herb.mat, F))
beta.dpw$hd0.null <- rowMeans(null)
beta.dpw$hd0.ses <- (beta.dpw$hd0 - beta.dpw$hd0.null)/apply(null, 1, sd)

beta.dpw$hd1 <- as.vector(comdist(ha, herb.mat, T))
null <- replicate(100, dpw.rand(ha, herb.mat, T))
beta.dpw$hd1.null <- rowMeans(null)
beta.dpw$hd1.ses <- (beta.dpw$hd1 - beta.dpw$hd1.null)/apply(null, 1, sd)


beta.dpw$dh <- dist(hgt)
beta.dpw$mh <- as.dist(sapply(hgt, function(x) sapply(hgt, function(z) (x+z)/2)))

save(beta.dpw, file = "clean.data/beta.dpw.full.rda")