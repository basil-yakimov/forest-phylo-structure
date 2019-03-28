library(vegan)
library(picante)

load("clean.data/wood-phylo.rda")
load("clean.data/herb-phylo.rda")

wood.mat <- cophenetic(wood.tree)
herb.mat <- cophenetic(herb.tree)

load("clean.data/tree-abund.rda")
load("clean.data/shrub-abund.rda")
load("clean.data/herb-abund.rda")

hgt <- read.table("raw.data/hgt.txt")[[1]]

#---#

library(ecomf)

tpd <- PqD(ta, wood.tree, q = -1)

res <- oecosimu(round(ha), PqD, method = "quasiswap_count", nsimul = 999, burnin = 100, thin = 100,
                tree = herb.tree, q = 2)

plot(hgt, res$oecosimu$z)


hist(res$oecosimu$simulated[95, ])


tar <- randomizeMatrix(ta, null.model = "independentswap")
cbind(colSums(ta), colSums(tar))
cbind(rowSums(ta), rowSums(tar))
cbind(colSums(ta > 0), colSums(tar > 0))
cbind(rowSums(ta > 0), rowSums(tar > 0))


null <- replicate(1000, mpd(randomizeMatrix(ta, null.model = "independentswap"), wood.mat))
ta.ses.mpd <- (mpd(ta, wood.mat) - rowMeans(null)) / apply(null, 1, sd)

null <- replicate(1000, mpd(randomizeMatrix(ta, null.model = "independentswap"), wood.mat, abundance.weighted = T))
ta.ses.mpda <- (mpd(ta, wood.mat, abundance.weighted = T) - rowMeans(null)) / apply(null, 1, sd)

null <- replicate(1000, PqD(randomizeMatrix(ta, null.model = "independentswap"), wood.tree, 0, F))
ta.ses.pd0 <- (PqD(ta, wood.tree, q = 0, hill = F) - rowMeans(null)) / apply(null, 1, sd)

null <- replicate(1000, PqD(randomizeMatrix(ta, null.model = "independentswap"), wood.tree, 1, F))
ta.ses.pd1 <- (PqD(ta, wood.tree, q = 1, hill = F) - rowMeans(null)) / apply(null, 1, sd)


op <- par(mfrow = c(2, 2))
plot(hgt, ta.ses.mpd, pch = 19, col = "tomato", main = "MPD")
plot(hgt, ta.ses.mpda, pch = 19, col = "tomato", main = "MPDa")
plot(hgt, ta.ses.pd0, pch = 19, col = "tomato", main = expression(PD[0]))
plot(hgt, ta.ses.pd1, pch = 19, col = "tomato", main = expression(PD[1]))
par(op)


null <- replicate(1000, mpd(randomizeMatrix(sa, null.model = "independentswap"), wood.mat))
sa.ses.mpd <- (mpd(sa, wood.mat) - rowMeans(null)) / apply(null, 1, sd)

null <- replicate(1000, mpd(randomizeMatrix(sa, null.model = "independentswap"), wood.mat, abundance.weighted = T))
sa.ses.mpda <- (mpd(sa, wood.mat, abundance.weighted = T) - rowMeans(null)) / apply(null, 1, sd)

null <- replicate(1000, PqD(randomizeMatrix(sa, null.model = "independentswap"), wood.tree, 0, F))
sa.ses.pd0 <- (PqD(sa, wood.tree, q = 0, hill = F) - rowMeans(null)) / apply(null, 1, sd)

null <- replicate(1000, PqD(randomizeMatrix(sa, null.model = "independentswap"), wood.tree, 1, F))
sa.ses.pd1 <- (PqD(sa, wood.tree, q = 1, hill = F) - rowMeans(null)) / apply(null, 1, sd)


op <- par(mfrow = c(2, 2))
plot(hgt, sa.ses.mpd, pch = 19, col = "skyblue", main = "MPD")
plot(hgt, sa.ses.mpda, pch = 19, col = "skyblue", main = "MPDa")
plot(hgt, sa.ses.pd0, pch = 19, col = "skyblue", main = expression(PD[0]))
plot(hgt, sa.ses.pd1, pch = 19, col = "skyblue", main = expression(PD[1]))
par(op)


null <- replicate(1000, mpd(randomizeMatrix(ha, null.model = "independentswap"), herb.mat))
ha.ses.mpd <- (mpd(ha, herb.mat) - rowMeans(null)) / apply(null, 1, sd)

null <- replicate(1000, mpd(randomizeMatrix(ha, null.model = "independentswap"), herb.mat, abundance.weighted = T))
ha.ses.mpda <- (mpd(ha, herb.mat, abundance.weighted = T) - rowMeans(null)) / apply(null, 1, sd)

null <- replicate(1000, PqD(randomizeMatrix(ha, null.model = "independentswap"), herb.tree, 0, F))
ha.ses.pd0 <- (PqD(ha, herb.tree, q = 0, hill = F) - rowMeans(null)) / apply(null, 1, sd)

null <- replicate(1000, PqD(randomizeMatrix(ha, null.model = "independentswap"), herb.tree, 1, F))
ha.ses.pd1 <- (PqD(ha, herb.tree, q = 1, hill = F) - rowMeans(null)) / apply(null, 1, sd)


op <- par(mfrow = c(2, 2))
plot(hgt, ha.ses.mpd, pch = 19, col = "forestgreen", main = "MPD")
plot(hgt, ha.ses.mpda, pch = 19, col = "forestgreen", main = "MPDa")
plot(hgt, ha.ses.pd0, pch = 19, col = "forestgreen", main = expression(PD[0]))
plot(hgt, ha.ses.pd1, pch = 19, col = "forestgreen", main = expression(PD[1]))
par(op)

null <- replicate(1000, PqD(randomizeMatrix(ta, null.model = "independentswap"), wood.tree, 2, F))
ta.ses.pd2 <- (PqD(ta, wood.tree, q = 2, hill = F) - rowMeans(null)) / apply(null, 1, sd)

null <- replicate(1000, PqD(randomizeMatrix(sa, null.model = "independentswap"), wood.tree, 2, F))
sa.ses.pd2 <- (PqD(sa, wood.tree, q = 2, hill = F) - rowMeans(null)) / apply(null, 1, sd)

null <- replicate(1000, PqD(randomizeMatrix(ha, null.model = "independentswap"), herb.tree, 2, F))
ha.ses.pd2 <- (PqD(ha, herb.tree, q = 2, hill = F) - rowMeans(null)) / apply(null, 1, sd)


null <- replicate(1000, PqD(randomizeMatrix(ta, null.model = "independentswap"), wood.tree, -1, F))
ta.ses.pdn1 <- (PqD(ta, wood.tree, q = -1, hill = F) - rowMeans(null)) / apply(null, 1, sd)

null <- replicate(1000, PqD(randomizeMatrix(sa, null.model = "independentswap"), wood.tree, -1, F))
sa.ses.pdn1 <- (PqD(sa, wood.tree, q = -1, hill = F) - rowMeans(null)) / apply(null, 1, sd)

null <- replicate(1000, PqD(randomizeMatrix(ha, null.model = "independentswap"), herb.tree, -1, F))
ha.ses.pdn1 <- (PqD(ha, herb.tree, q = -1, hill = F) - rowMeans(null)) / apply(null, 1, sd)


plot.ses <- function(ses, col, lab) {
  plot(hgt, ses, pch = 19, col = col, ylab = lab)
  fit <- lm(ses ~ hgt)
  if (anova(fit)[1, 5] < 0.05) {
    abline(fit)
    usr <- par("usr")
    text(usr[1], usr[3] + (usr[4] - usr[3]) / 10, round(cor(ses, hgt, use = "complete")^2, digits = 3), pos = 4, adj = c(1, 1))
  }
}


op <- par(mfcol = c(4, 3))

plot.ses(ta.ses.pdn1, col = "tomato", lab = expression(PD[-1]))
plot.ses(ta.ses.pd0, col = "tomato", lab = expression(PD[0]))
plot.ses(ta.ses.pd1, col = "tomato", lab = expression(PD[1]))
plot.ses(ta.ses.pd2, col = "tomato", lab = expression(PD[2]))


plot.ses(sa.ses.pdn1, col = "skyblue", lab = expression(PD[-1]))
plot.ses(sa.ses.pd0, col = "skyblue", lab = expression(PD[0]))
plot.ses(sa.ses.pd1, col = "skyblue", lab = expression(PD[1]))
plot.ses(sa.ses.pd2, col = "skyblue", lab = expression(PD[2]))

plot.ses(ha.ses.pdn1, col = "forestgreen", lab = expression(PD[-1]))
plot.ses(ha.ses.pd0, col = "forestgreen", lab = expression(PD[0]))
plot.ses(ha.ses.pd1, col = "forestgreen", lab = expression(PD[1]))
plot.ses(ha.ses.pd2, col = "forestgreen", lab = expression(PD[2]))

par(op)