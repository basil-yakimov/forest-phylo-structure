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

#---tree layer analysis---#

null <- replicate(999, mpd(randomizeMatrix(ta, null.model = "independentswap"), wood.mat))
ta.ses <- data.frame(mpd = (mpd(ta, wood.mat) - rowMeans(null)) / apply(null, 1, sd))

null <- replicate(999, mpd(randomizeMatrix(ta, null.model = "independentswap"), wood.mat, ab = T))
ta.ses$mpd.a <- (mpd(ta, wood.mat, ab = T) - rowMeans(null)) / apply(null, 1, sd)

null <- replicate(1000, PqD(randomizeMatrix(ta, null.model = "independentswap"), wood.tree, -1, F))
ta.ses$pqd.n1 <- (PqD(ta, wood.tree, q = -1, hill = F) - rowMeans(null)) / apply(null, 1, sd)

null <- replicate(1000, PqD(randomizeMatrix(ta, null.model = "independentswap"), wood.tree, 0, F))
ta.ses$pqd.0 <- (PqD(ta, wood.tree, q = 0, hill = F) - rowMeans(null)) / apply(null, 1, sd)

null <- replicate(1000, PqD(randomizeMatrix(ta, null.model = "independentswap"), wood.tree, 1, F))
ta.ses$pqd.1 <- (PqD(ta, wood.tree, q = 1, hill = F) - rowMeans(null)) / apply(null, 1, sd)

null <- replicate(1000, PqD(randomizeMatrix(ta, null.model = "independentswap"), wood.tree, 2, F))
ta.ses$pqd.2 <- (PqD(ta, wood.tree, q = 2, hill = F) - rowMeans(null)) / apply(null, 1, sd)

null <- replicate(1000, PqD(randomizeMatrix(ta, null.model = "independentswap"), wood.tree, 5, F))
ta.ses$pqd.5 <- (PqD(ta, wood.tree, q = 5, hill = F) - rowMeans(null)) / apply(null, 1, sd)

#---shrub layer analysis---#

null <- replicate(999, mpd(randomizeMatrix(sa, null.model = "independentswap"), shrub.mat))
sa.ses <- data.frame(mpd = (mpd(sa, shrub.mat) - rowMeans(null)) / apply(null, 1, sd))

null <- replicate(999, mpd(randomizeMatrix(sa, null.model = "independentswap"), shrub.mat, ab = T))
sa.ses$mpd.a <- (mpd(sa, shrub.mat, ab = T) - rowMeans(null)) / apply(null, 1, sd)

null <- replicate(1000, PqD(randomizeMatrix(sa, null.model = "independentswap"), shrub.tree, -1, F))
sa.ses$pqd.n1 <- (PqD(sa, shrub.tree, q = -1, hill = F) - rowMeans(null)) / apply(null, 1, sd)

null <- replicate(1000, PqD(randomizeMatrix(sa, null.model = "independentswap"), shrub.tree, 0, F))
sa.ses$pqd.0 <- (PqD(sa, shrub.tree, q = 0, hill = F) - rowMeans(null)) / apply(null, 1, sd)

null <- replicate(1000, PqD(randomizeMatrix(sa, null.model = "independentswap"), shrub.tree, 1, F))
sa.ses$pqd.1 <- (PqD(sa, shrub.tree, q = 1, hill = F) - rowMeans(null)) / apply(null, 1, sd)

null <- replicate(1000, PqD(randomizeMatrix(sa, null.model = "independentswap"), shrub.tree, 2, F))
sa.ses$pqd.2 <- (PqD(sa, shrub.tree, q = 2, hill = F) - rowMeans(null)) / apply(null, 1, sd)

null <- replicate(1000, PqD(randomizeMatrix(sa, null.model = "independentswap"), shrub.tree, 5, F))
sa.ses$pqd.5 <- (PqD(sa, shrub.tree, q = 5, hill = F) - rowMeans(null)) / apply(null, 1, sd)

#---herb layer analysis---#

null <- replicate(999, mpd(randomizeMatrix(ha, null.model = "independentswap"), herb.mat))
ha.ses <- data.frame(mpd = (mpd(ha, herb.mat) - rowMeans(null)) / apply(null, 1, sd))

null <- replicate(999, mpd(randomizeMatrix(ha, null.model = "independentswap"), herb.mat, ab = T))
ha.ses$mpd.a <- (mpd(ha, herb.mat, ab = T) - rowMeans(null)) / apply(null, 1, sd)

null <- replicate(1000, PqD(randomizeMatrix(ha, null.model = "independentswap"), herb.tree, -1, F))
ha.ses$pqd.n1 <- (PqD(ha, herb.tree, q = -1, hill = F) - rowMeans(null)) / apply(null, 1, sd)

null <- replicate(1000, PqD(randomizeMatrix(ha, null.model = "independentswap"), herb.tree, 0, F))
ha.ses$pqd.0 <- (PqD(ha, herb.tree, q = 0, hill = F) - rowMeans(null)) / apply(null, 1, sd)

null <- replicate(1000, PqD(randomizeMatrix(ha, null.model = "independentswap"), herb.tree, 1, F))
ha.ses$pqd.1 <- (PqD(ha, herb.tree, q = 1, hill = F) - rowMeans(null)) / apply(null, 1, sd)

null <- replicate(1000, PqD(randomizeMatrix(ha, null.model = "independentswap"), herb.tree, 2, F))
ha.ses$pqd.2 <- (PqD(ha, herb.tree, q = 2, hill = F) - rowMeans(null)) / apply(null, 1, sd)

null <- replicate(1000, PqD(randomizeMatrix(ha, null.model = "independentswap"), herb.tree, 5, F))
ha.ses$pqd.5 <- (PqD(ha, herb.tree, q = 5, hill = F) - rowMeans(null)) / apply(null, 1, sd)

#--- add hgt and save results ---#

hgt <- read.table("raw.data/hgt.txt")[[1]]
ta.ses$hgt <- sa.ses$hgt <- ha.ses$hgt <- hgt

save(ta.ses, sa.ses, ha.ses, file = "clean.data/phylo-ses-pqd.rda")
