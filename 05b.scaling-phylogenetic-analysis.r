library(ape)
library(picante)
library(ecomf)

load("clean.data/tree-abund.rda")
load("clean.data/shrub-abund.rda")
load("clean.data/herb-abund.rda")

load("clean.data/wood-phylo.rda")
load("clean.data/herb-phylo.rda")

shrub.tree <- drop.tip(wood.tree, wood.tree$tip.label[!(wood.tree$tip.label %in% colnames(sa))])
wood.tree <- drop.tip(wood.tree, wood.tree$tip.label[!(wood.tree$tip.label %in% colnames(ta))])

wood.mat <- cophenetic(wood.tree)
shrub.mat <- cophenetic(shrub.tree)
herb.mat <- cophenetic(herb.tree)

load("clean.data/scaling.rda")

#________________________________________________

ta.ses.sc <- sa.ses.sc <- ha.ses.sc <- vector(mode = "list", length = 20)

null <- "taxa.labels"

for (jj in 1:20)
{
  print(jj)

  # tree
  
  null <- replicate(999, mpd(randomizeMatrix(ta_sc[[jj]], null.model = "independentswap"), wood.mat))
  ta.ses.sc[[jj]] <- data.frame(mpd = (mpd(ta_sc[[jj]], wood.mat) - rowMeans(null)) / apply(null, 1, sd))
  
  null <- replicate(999, mpd(randomizeMatrix(ta_sc[[jj]], null.model = "independentswap"), wood.mat, ab = T))
  ta.ses.sc[[jj]]$mpd.a <- (mpd(ta_sc[[jj]], wood.mat, ab = T) - rowMeans(null)) / apply(null, 1, sd)
  
  null <- replicate(1000, PqD(randomizeMatrix(ta_sc[[jj]], null.model = "independentswap"), wood.tree, -1, F))
  ta.ses.sc[[jj]]$pqd.n1 <- (PqD(ta_sc[[jj]], wood.tree, q = -1, hill = F) - rowMeans(null)) / apply(null, 1, sd)
  
  null <- replicate(1000, PqD(randomizeMatrix(ta_sc[[jj]], null.model = "independentswap"), wood.tree, 0, F))
  ta.ses.sc[[jj]]$pqd.0 <- (PqD(ta_sc[[jj]], wood.tree, q = 0, hill = F) - rowMeans(null)) / apply(null, 1, sd)
  
  null <- replicate(1000, PqD(randomizeMatrix(ta_sc[[jj]], null.model = "independentswap"), wood.tree, 1, F))
  ta.ses.sc[[jj]]$pqd.1 <- (PqD(ta_sc[[jj]], wood.tree, q = 1, hill = F) - rowMeans(null)) / apply(null, 1, sd)
  
  null <- replicate(1000, PqD(randomizeMatrix(ta_sc[[jj]], null.model = "independentswap"), wood.tree, 2, F))
  ta.ses.sc[[jj]]$pqd.2 <- (PqD(ta_sc[[jj]], wood.tree, q = 2, hill = F) - rowMeans(null)) / apply(null, 1, sd)
  
  null <- replicate(1000, PqD(randomizeMatrix(ta_sc[[jj]], null.model = "independentswap"), wood.tree, 5, F))
  ta.ses.sc[[jj]]$pqd.5 <- (PqD(ta_sc[[jj]], wood.tree, q = 5, hill = F) - rowMeans(null)) / apply(null, 1, sd)
  
  # shrub
  
  null <- replicate(999, mpd(randomizeMatrix(sa_sc[[jj]], null.model = "independentswap"), shrub.mat))
  sa.ses.sc[[jj]] <- data.frame(mpd = (mpd(sa_sc[[jj]], shrub.mat) - rowMeans(null)) / apply(null, 1, sd))
  
  null <- replicate(999, mpd(randomizeMatrix(sa_sc[[jj]], null.model = "independentswap"), shrub.mat, ab = T))
  sa.ses.sc[[jj]]$mpd.a <- (mpd(sa_sc[[jj]], shrub.mat, ab = T) - rowMeans(null)) / apply(null, 1, sd)
  
  null <- replicate(1000, PqD(randomizeMatrix(sa_sc[[jj]], null.model = "independentswap"), shrub.tree, -1, F))
  sa.ses.sc[[jj]]$pqd.n1 <- (PqD(sa_sc[[jj]], shrub.tree, q = -1, hill = F) - rowMeans(null)) / apply(null, 1, sd)
  
  null <- replicate(1000, PqD(randomizeMatrix(sa_sc[[jj]], null.model = "independentswap"), shrub.tree, 0, F))
  sa.ses.sc[[jj]]$pqd.0 <- (PqD(sa_sc[[jj]], shrub.tree, q = 0, hill = F) - rowMeans(null)) / apply(null, 1, sd)
  
  null <- replicate(1000, PqD(randomizeMatrix(sa_sc[[jj]], null.model = "independentswap"), shrub.tree, 1, F))
  sa.ses.sc[[jj]]$pqd.1 <- (PqD(sa_sc[[jj]], shrub.tree, q = 1, hill = F) - rowMeans(null)) / apply(null, 1, sd)
  
  null <- replicate(1000, PqD(randomizeMatrix(sa_sc[[jj]], null.model = "independentswap"), shrub.tree, 2, F))
  sa.ses.sc[[jj]]$pqd.2 <- (PqD(sa_sc[[jj]], shrub.tree, q = 2, hill = F) - rowMeans(null)) / apply(null, 1, sd)
  
  null <- replicate(1000, PqD(randomizeMatrix(sa_sc[[jj]], null.model = "independentswap"), shrub.tree, 5, F))
  sa.ses.sc[[jj]]$pqd.5 <- (PqD(sa_sc[[jj]], shrub.tree, q = 5, hill = F) - rowMeans(null)) / apply(null, 1, sd)
  
  # herb
  
  null <- replicate(999, mpd(randomizeMatrix(ha_sc[[jj]], null.model = "independentswap"), herb.mat))
  ha.ses.sc[[jj]] <- data.frame(mpd = (mpd(ha_sc[[jj]], herb.mat) - rowMeans(null)) / apply(null, 1, sd))
  
  null <- replicate(999, mpd(randomizeMatrix(ha_sc[[jj]], null.model = "independentswap"), herb.mat, ab = T))
  ha.ses.sc[[jj]]$mpd.a <- (mpd(ha_sc[[jj]], herb.mat, ab = T) - rowMeans(null)) / apply(null, 1, sd)
  
  null <- replicate(1000, PqD(randomizeMatrix(ha_sc[[jj]], null.model = "independentswap"), herb.tree, -1, F))
  ha.ses.sc[[jj]]$pqd.n1 <- (PqD(ha_sc[[jj]], herb.tree, q = -1, hill = F) - rowMeans(null)) / apply(null, 1, sd)
  
  null <- replicate(1000, PqD(randomizeMatrix(ha_sc[[jj]], null.model = "independentswap"), herb.tree, 0, F))
  ha.ses.sc[[jj]]$pqd.0 <- (PqD(ha_sc[[jj]], herb.tree, q = 0, hill = F) - rowMeans(null)) / apply(null, 1, sd)
  
  null <- replicate(1000, PqD(randomizeMatrix(ha_sc[[jj]], null.model = "independentswap"), herb.tree, 1, F))
  ha.ses.sc[[jj]]$pqd.1 <- (PqD(ha_sc[[jj]], herb.tree, q = 1, hill = F) - rowMeans(null)) / apply(null, 1, sd)
  
  null <- replicate(1000, PqD(randomizeMatrix(ha_sc[[jj]], null.model = "independentswap"), herb.tree, 2, F))
  ha.ses.sc[[jj]]$pqd.2 <- (PqD(ha_sc[[jj]], herb.tree, q = 2, hill = F) - rowMeans(null)) / apply(null, 1, sd)
  
  null <- replicate(1000, PqD(randomizeMatrix(ha_sc[[jj]], null.model = "independentswap"), herb.tree, 5, F))
  ha.ses.sc[[jj]]$pqd.5 <- (PqD(ha_sc[[jj]], herb.tree, q = 5, hill = F) - rowMeans(null)) / apply(null, 1, sd)
  
}  

save(ta.ses.sc, sa.ses.sc, ha.ses.sc, file = "clean.data/scaling-ses-pqd.rda")
