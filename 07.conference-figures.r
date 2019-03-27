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

png("figures/conf-tax-div.png", width = 3000, height = 2000)
op <- par(mfcol = c(3, 3), mar = c(4, 4, 0.5, 0.5), cex = 2)

plot(hgt, specnumber(ta), pch = 21, bg = "tomato", xlab = "Altitude, m", ylab = expression(D[0] == S))
legend("top", legend = "Tree layer", bty = "n")
plot(hgt, exp(diversity(ta)), pch = 21, bg = "tomato", xlab = "Altitude, m", ylab = expression(D[1]))
plot(hgt, diversity(ta, "invsimpson"), pch = 21, bg = "tomato", xlab = "Altitude, m", ylab = expression(D[2]))

plot(hgt, specnumber(sa), pch = 21, bg = "skyblue", xlab = "Altitude, m", ylab = expression(D[0] == S))
legend("top", legend = "Shrub layer", bty = "n")
plot(hgt, exp(diversity(sa)), pch = 21, bg = "skyblue", xlab = "Altitude, m", ylab = expression(D[1]))
plot(hgt, diversity(sa, "invsimpson"), pch = 21, bg = "skyblue", xlab = "Altitude, m", ylab = expression(D[2]))

plot(hgt, specnumber(ha), pch = 21, bg = "forestgreen", xlab = "Altitude, m", ylab = expression(D[0] == S))
legend("top", legend = "Herb layer", bty = "n")
plot(hgt, exp(diversity(ha)), pch = 21, bg = "forestgreen", xlab = "Altitude, m", ylab = expression(D[1]))
plot(hgt, diversity(ha, "invsimpson"), pch = 21, bg = "forestgreen", xlab = "Altitude, m", ylab = expression(D[2]))

par(op)
dev.off()

#---#

library(ecomf)

tpd <- PqD(ta, wood.tree)
spd <- PqD(sa, wood.tree)
hpd <- PqD(ha, herb.tree)

png("figures/conf-phylo-div.png", width = 3000, height = 2000)
op <- par(mfcol = c(3, 3), mar = c(4, 4, 0.5, 0.5), cex = 2)

plot(hgt, tpd[, 1], pch = 21, bg = "tomato", xlab = "Altitude, m", ylab = expression(PD[0]))
legend("top", legend = "Tree layer", bty = "n")
plot(hgt, tpd[, 2], pch = 21, bg = "tomato", xlab = "Altitude, m", ylab = expression(PD[1]))
plot(hgt, tpd[, 3], pch = 21, bg = "tomato", xlab = "Altitude, m", ylab = expression(PD[2]))

plot(hgt, spd[, 1], pch = 21, bg = "skyblue", xlab = "Altitude, m", ylab = expression(PD[0]))
legend("top", legend = "Shrub layer", bty = "n")
plot(hgt, spd[, 2], pch = 21, bg = "skyblue", xlab = "Altitude, m", ylab = expression(PD[1]))
plot(hgt, spd[, 3], pch = 21, bg = "skyblue", xlab = "Altitude, m", ylab = expression(PD[2]))

plot(hgt, hpd[, 1], pch = 21, bg = "forestgreen", xlab = "Altitude, m", ylab = expression(PD[0]))
legend("top", legend = "Herb layer", bty = "n")
plot(hgt, hpd[, 2], pch = 21, bg = "forestgreen", xlab = "Altitude, m", ylab = expression(PD[1]))
plot(hgt, hpd[, 3], pch = 21, bg = "forestgreen", xlab = "Altitude, m", ylab = expression(PD[2]))

par(op)
dev.off()