library(vegan)

source("R/plot.ses.r")

load("clean.data/wood-phylo.rda")
load("clean.data/herb-phylo.rda")

wood.mat <- cophenetic(wood.tree)
herb.mat <- cophenetic(herb.tree)

load("clean.data/tree-abund.rda")
load("clean.data/shrub-abund.rda")
load("clean.data/herb-abund.rda")

hgt <- read.table("raw.data/hgt.txt")[[1]]

#---#

png("figures/semifinal-Fig01-tax-div.png", width = 3000, height = 1500)
op <- par(mfcol = c(2, 3), mar = c(4, 4, 0.5, 0.5), cex = 2)

plot.ses(specnumber(ta), hgt, col = "tomato", lab = expression(D[0] == S), xlab = "Altitude, m")
legend("top", legend = "Tree layer", bty = "n")
plot.ses(diversity(ta, "invsimpson"), hgt, col = "tomato", lab = expression(D[2]), xlab = "Altitude, m")

plot.ses(specnumber(sa), hgt, col = "skyblue", lab = expression(D[0] == S), xlab = "Altitude, m")
legend("top", legend = "Shrub layer", bty = "n")
plot.ses(diversity(sa, "invsimpson"), hgt, col = "skyblue", lab = expression(D[2]), xlab = "Altitude, m")

plot.ses(specnumber(ha), hgt, col = "forestgreen", lab = expression(D[0] == S), xlab = "Altitude, m")
legend("top", legend = "Herb layer", bty = "n")
plot.ses(diversity(ha, "invsimpson"), hgt, col = "forestgreen", lab = expression(D[2]), xlab = "Altitude, m")

par(op)
dev.off()

#---#

library(ecomf)

tpd <- PqD(ta, wood.tree)
spd <- PqD(sa, wood.tree)
hpd <- PqD(ha, herb.tree)

png("figures/semifinal-Fig02-phylo-div.png", width = 3000, height = 1500)
op <- par(mfcol = c(2, 3), mar = c(4, 4, 0.5, 0.5), cex = 2)

plot.ses(tpd[, 1], hgt, col = "tomato", lab = expression(PD[0]), xlab = "Altitude, m")
legend("top", legend = "Tree layer", bty = "n")
plot.ses(tpd[, 3], hgt, col = "tomato", lab = expression(PD[2]), xlab = "Altitude, m")

plot.ses(spd[, 1], hgt, col = "skyblue", lab = expression(PD[0]), xlab = "Altitude, m")
legend("top", legend = "Shrub layer", bty = "n")
plot.ses(spd[, 3], hgt, col = "skyblue", lab = expression(PD[2]), xlab = "Altitude, m")

plot.ses(hpd[, 1], hgt, col = "forestgreen", lab = expression(PD[0]), xlab = "Altitude, m")
legend("top", legend = "Herb layer", bty = "n")
plot.ses(hpd[, 3], hgt, col = "forestgreen", lab = expression(PD[2]), xlab = "Altitude, m")

par(op)
dev.off()

#---#

load("clean.data/phylo-ses-pqd.rda")

png("figures/semifinal-Fig03-NRI.png", width = 3000, height = 1500)
op <- par(mfcol = c(2, 3), mar = c(4, 4, 0.5, 0.5), cex = 2)

plot.ses(ta.ses$mpd, hgt, col = "tomato", lab = "-NRI", xlab = "Altitude, m")
legend("topright", legend = "Tree layer", bty = "n")
plot.ses(ta.ses$mpd.a, hgt, col = "tomato", lab = expression(-NRI[a]), xlab = "Altitude, m")

plot.ses(sa.ses$mpd, hgt, col = "skyblue", lab = "-NRI", xlab = "Altitude, m")
legend("top", legend = "Shrub layer", bty = "n")
plot.ses(sa.ses$mpd.a, hgt, col = "skyblue", lab = expression(-NRI[a]), xlab = "Altitude, m")

plot.ses(ha.ses$mpd, hgt, col = "forestgreen", lab = "-NRI", xlab = "Altitude, m")
legend("top", legend = "Herb layer", bty = "n")
plot.ses(ha.ses$mpd.a, hgt, col = "forestgreen", lab = expression(-NRI[a]), xlab = "Altitude, m")

par(op)
dev.off()


png("figures/semifinal-Fig04-NTI.png", width = 3000, height = 1500)
op <- par(mfcol = c(2, 3), mar = c(4, 4, 0.5, 0.5), cex = 2)

plot.ses(ta.ses$mntd, hgt, col = "tomato", lab = "-NTI", xlab = "Altitude, m")
legend("topright", legend = "Tree layer", bty = "n")
plot.ses(ta.ses$mntd.a, hgt, col = "tomato", lab = expression(-NTI[a]), xlab = "Altitude, m")

plot.ses(sa.ses$mntd, hgt, col = "skyblue", lab = "-NTI", xlab = "Altitude, m")
legend("top", legend = "Shrub layer", bty = "n")
plot.ses(sa.ses$mntd.a, hgt, col = "skyblue", lab = expression(-NTI[a]), xlab = "Altitude, m")

plot.ses(ha.ses$mntd, hgt, col = "forestgreen", lab = "-NTI", xlab = "Altitude, m")
legend("top", legend = "Herb layer", bty = "n")
plot.ses(ha.ses$mntd.a, hgt, col = "forestgreen", lab = expression(-NTI[a]), xlab = "Altitude, m")

par(op)
dev.off()

#---#

dh <- dist(hgt)

ta.sor <- vegdist(ta, method = "bray", binary = T)
ta.mh <- vegdist(ta, method = "horn", binary = F)
sa.sor <- vegdist(sa, method = "bray", binary = T)
sa.mh <- vegdist(sa, method = "horn", binary = F)
ha.sor <- vegdist(ha, method = "bray", binary = T)
ha.mh <- vegdist(ha, method = "horn", binary = F)

png("figures/semifinal-Fig05-tax-beta-dh.png", width = 3000, height = 1500)
op <- par(mfcol = c(2, 3), mar = c(4, 4, 0.5, 0.5), cex = 2)

plot.ses(ta.sor, dh, col = "tomato", lab = "Sorensen index", xlab = "Altitude difference, m")
legend("topright", legend = "Tree layer", bty = "n")
plot.ses(ta.mh, dh, col = "tomato", lab = "Morisita-Horn index", xlab = "Altitude difference, m")

plot.ses(sa.sor, dh, col = "skyblue", lab = "Sorensen index", xlab = "Altitude difference, m")
legend("topright", legend = "Shrub layer", bty = "n")
plot.ses(sa.mh, dh, col = "skyblue", lab = "Morisita-Horn index", xlab = "Altitude difference, m")

plot.ses(ha.sor, dh, col = "forestgreen", lab = "Sorensen index", xlab = "Altitude difference, m")
legend("topright", legend = "Herb layer", bty = "n")
plot.ses(ha.mh, dh, col = "forestgreen", lab = "Morisita-Horn index", xlab = "Altitude difference, m")

par(op)
dev.off()

ta.sor.95 <- diag(as.matrix(ta.sor)[-1, ])
ta.mh.95 <- diag(as.matrix(ta.mh)[-1, ])
sa.sor.95 <- diag(as.matrix(sa.sor)[-1, ])
sa.mh.95 <- diag(as.matrix(sa.mh)[-1, ])
ha.sor.95 <- diag(as.matrix(ha.sor)[-1, ])
ha.mh.95 <- diag(as.matrix(ha.mh)[-1, ])

mh <- (hgt[1:95] + hgt[2:96])/2

png("figures/semifinal-Fig06-tax-beta-mh.png", width = 3000, height = 1500)
op <- par(mfcol = c(2, 3), mar = c(4, 4, 0.5, 0.5), cex = 2)

plot.ses(ta.sor.95[c(T, F)], mh[c(T, F)], col = "tomato", lab = "Sorensen index", xlab = "Mean altitude, m")
legend("topright", legend = "Tree layer", bty = "n")
plot.ses(ta.mh.95[c(T, F)], mh[c(T, F)], col = "tomato", lab = "Morisita-Horn index", xlab = "Mean altitude, m")

plot.ses(sa.sor.95[c(T, F)], mh[c(T, F)], col = "skyblue", lab = "Sorensen index", xlab = "Mean altitude, m")
legend("topright", legend = "Shrub layer", bty = "n")
plot.ses(sa.mh.95[c(T, F)], mh[c(T, F)], col = "skyblue", lab = "Morisita-Horn index", xlab = "Mean altitude, m")

plot.ses(ha.sor.95[c(T, F)], mh[c(T, F)], col = "forestgreen", lab = "Sorensen index", xlab = "Mean altitude, m")
legend("topright", legend = "Herb layer", bty = "n")
plot.ses(ha.mh.95[c(T, F)], mh[c(T, F)], col = "forestgreen", lab = "Morisita-Horn index", xlab = "Mean altitude, m")

par(op)
dev.off()

#---#

library(picante)

ta.dpw0 <- comdist(ta, wood.mat, F)
ta.dpw1 <- comdist(ta, wood.mat, T)
sa.dpw0 <- comdist(sa, wood.mat, F)
sa.dpw1 <- comdist(sa, wood.mat, T)
ha.dpw0 <- comdist(ha, herb.mat, F)
ha.dpw1 <- comdist(ha, herb.mat, T)

ta.dnn0 <- comdistnt(ta, wood.mat, F)
ta.dnn1 <- comdistnt(ta, wood.mat, T)
sa.dnn0 <- comdistnt(sa, wood.mat, F)
sa.dnn1 <- comdistnt(sa, wood.mat, T)
ha.dnn0 <- comdistnt(ha, herb.mat, F)
ha.dnn1 <- comdistnt(ha, herb.mat, T)

png("figures/semifinal-Fig07-dpw-dh.png", width = 3000, height = 1500)
op <- par(mfcol = c(2, 3), mar = c(4, 4, 0.5, 0.5), cex = 2)

plot.ses(ta.dpw0, dh, col = "tomato", lab = expression(D[pw]), xlab = "Altitude difference, m")
legend("topright", legend = "Tree layer", bty = "n")
plot.ses(ta.dpw1, dh, col = "tomato", lab = expression("D'"[pw]), xlab = "Altitude difference, m")

plot.ses(sa.dpw0, dh, col = "skyblue", lab = expression(D[pw]), xlab = "Altitude difference, m")
legend("topright", legend = "Shrub layer", bty = "n")
plot.ses(sa.dpw1, dh, col = "skyblue", lab = expression("D'"[pw]), xlab = "Altitude difference, m")

plot.ses(ha.dpw0, dh, col = "forestgreen", lab = expression(D[pw]), xlab = "Altitude difference, m")
legend("topright", legend = "Herb layer", bty = "n")
plot.ses(ha.dpw1, dh, col = "forestgreen", lab = expression("D'"[pw]), xlab = "Altitude difference, m")

par(op)
dev.off()

ta.dpw0.95 <- diag(as.matrix(ta.dpw0)[-1, ])
ta.dpw1.95 <- diag(as.matrix(ta.dpw1)[-1, ])
sa.dpw0.95 <- diag(as.matrix(sa.dpw0)[-1, ])
sa.dpw1.95 <- diag(as.matrix(sa.dpw1)[-1, ])
ha.dpw0.95 <- diag(as.matrix(ha.dpw0)[-1, ])
ha.dpw1.95 <- diag(as.matrix(ha.dpw1)[-1, ])

ta.dnn0.95 <- diag(as.matrix(ta.dnn0)[-1, ])
ta.dnn1.95 <- diag(as.matrix(ta.dnn1)[-1, ])
sa.dnn0.95 <- diag(as.matrix(sa.dnn0)[-1, ])
sa.dnn1.95 <- diag(as.matrix(sa.dnn1)[-1, ])
ha.dnn0.95 <- diag(as.matrix(ha.dnn0)[-1, ])
ha.dnn1.95 <- diag(as.matrix(ha.dnn1)[-1, ])

png("figures/semifinal-Fig08-dpw-mh.png", width = 3000, height = 1500)
op <- par(mfcol = c(2, 3), mar = c(4, 4, 0.5, 0.5), cex = 2)

plot.ses(ta.dpw0.95[c(T, F)], mh[c(T, F)], col = "tomato", lab = expression(D[pw]), xlab = "Mean altitude, m")
legend("topright", legend = "Tree layer", bty = "n")
plot.ses(ta.dpw1.95[c(T, F)], mh[c(T, F)], col = "tomato", lab = expression("D'"[pw]), xlab = "Mean altitude, m")

plot.ses(sa.dpw0.95[c(T, F)], mh[c(T, F)], col = "skyblue", lab = expression(D[pw]), xlab = "Mean altitude, m")
legend("topright", legend = "Shrub layer", bty = "n")
plot.ses(sa.dpw1.95[c(T, F)], mh[c(T, F)], col = "skyblue", lab = expression("D'"[pw]), xlab = "Mean altitude, m")

plot.ses(ha.dpw0.95[c(T, F)], mh[c(T, F)], col = "forestgreen", lab = expression(D[pw]), xlab = "Mean altitude, m")
legend("topright", legend = "Herb layer", bty = "n")
plot.ses(ha.dpw1.95[c(T, F)], mh[c(T, F)], col = "forestgreen", lab = expression("D'"[pw]), xlab = "Mean altitude, m")

par(op)
dev.off()

png("figures/semifinal-Fig09-dnn-dh.png", width = 3000, height = 1500)
op <- par(mfcol = c(2, 3), mar = c(4, 4, 0.5, 0.5), cex = 2)

plot.ses(ta.dnn0, dh, col = "tomato", lab = expression(D[nn]), xlab = "Altitude difference, m")
legend("topright", legend = "Tree layer", bty = "n")
plot.ses(ta.dnn1, dh, col = "tomato", lab = expression("D'"[nn]), xlab = "Altitude difference, m")

plot.ses(sa.dnn0, dh, col = "skyblue", lab = expression(D[nn]), xlab = "Altitude difference, m")
legend("topright", legend = "Shrub layer", bty = "n")
plot.ses(sa.dnn1, dh, col = "skyblue", lab = expression("D'"[nn]), xlab = "Altitude difference, m")

plot.ses(ha.dnn0, dh, col = "forestgreen", lab = expression(D[nn]), xlab = "Altitude difference, m")
legend("topright", legend = "Herb layer", bty = "n")
plot.ses(ha.dnn1, dh, col = "forestgreen", lab = expression("D'"[nn]), xlab = "Altitude difference, m")

par(op)
dev.off()

png("figures/semifinal-Fig10-dnn-mh.png", width = 3000, height = 1500)
op <- par(mfcol = c(2, 3), mar = c(4, 4, 0.5, 0.5), cex = 2)

plot.ses(ta.dnn0.95[c(T, F)], mh[c(T, F)], col = "tomato", lab = expression(D[nn]), xlab = "Mean altitude, m")
legend("topright", legend = "Tree layer", bty = "n")
plot.ses(ta.dnn1.95[c(T, F)], mh[c(T, F)], col = "tomato", lab = expression("D'"[nn]), xlab = "Mean altitude, m")

plot.ses(sa.dnn0.95[c(T, F)], mh[c(T, F)], col = "skyblue", lab = expression(D[nn]), xlab = "Mean altitude, m")
legend("topright", legend = "Shrub layer", bty = "n")
plot.ses(sa.dnn1.95[c(T, F)], mh[c(T, F)], col = "skyblue", lab = expression("D'"[nn]), xlab = "Mean altitude, m")

plot.ses(ha.dnn0.95[c(T, F)], mh[c(T, F)], col = "forestgreen", lab = expression(D[nn]), xlab = "Mean altitude, m")
legend("topright", legend = "Herb layer", bty = "n")
plot.ses(ha.dnn1.95[c(T, F)], mh[c(T, F)], col = "forestgreen", lab = expression("D'"[nn]), xlab = "Mean altitude, m")

par(op)
dev.off()

#---#

load("clean.data/beta.dpw.full.rda")

png("figures/semifinal-Fig11-ses-dpw-dh.png", width = 3000, height = 1500)
op <- par(mfcol = c(2, 3), mar = c(4, 4, 0.5, 0.5), cex = 2)

plot.ses(beta.dpw$td0.ses, dh, col = "tomato", lab = expression("SES " * D[pw]), xlab = "Altitude difference, m")
legend("topright", legend = "Tree layer", bty = "n")
plot.ses(beta.dpw$td1.ses, dh, col = "tomato", lab = expression("SES D'"[pw]), xlab = "Altitude difference, m")

plot.ses(beta.dpw$sd0.ses, dh, col = "skyblue", lab = expression("SES " * D[pw]), xlab = "Altitude difference, m")
legend("topright", legend = "Shrub layer", bty = "n")
plot.ses(beta.dpw$sd1.ses, dh, col = "skyblue", lab = expression("SES D'"[pw]), xlab = "Altitude difference, m")

plot.ses(beta.dpw$hd0.ses, dh, col = "forestgreen", lab = expression("SES " * D[pw]), xlab = "Altitude difference, m")
legend("topright", legend = "Herb layer", bty = "n")
plot.ses(beta.dpw$hd1.ses, dh, col = "forestgreen", lab = expression("SES D'"[pw]), xlab = "Altitude difference, m")

par(op)
dev.off()


load("clean.data/pw-beta-tl.rda")


png("figures/semifinal-Fig12-ses-dpw-mh.png", width = 3000, height = 1500)
op <- par(mfcol = c(2, 3), mar = c(4, 4, 0.5, 0.5), cex = 2)

plot.ses(beta.ses$td0[c(T, F)], mh[c(T, F)], col = "tomato", lab = expression("SES " * D[pw]), xlab = "Mean altitude, m")
legend("topright", legend = "Tree layer", bty = "n")
plot.ses(beta.ses$td1[c(T, F)], mh[c(T, F)], col = "tomato", lab = expression("SES D'"[pw]), xlab = "Mean altitude, m")

plot.ses(beta.ses$sd0[c(T, F)], mh[c(T, F)], col = "skyblue", lab = expression("SES " * D[pw]), xlab = "Mean altitude, m")
legend("topright", legend = "Shrub layer", bty = "n")
plot.ses(beta.ses$sd1[c(T, F)], mh[c(T, F)], col = "skyblue", lab = expression("SES D'"[pw]), xlab = "Mean altitude, m")

plot.ses(beta.ses$hd0[c(T, F)], mh[c(T, F)], col = "forestgreen", lab = expression("SES " * D[pw]), xlab = "Mean altitude, m")
legend("topright", legend = "Herb layer", bty = "n")
plot.ses(beta.ses$hd1[c(T, F)], mh[c(T, F)], col = "forestgreen", lab = expression("SES D'"[pw]), xlab = "Mean altitude, m")

par(op)
dev.off()



png("figures/semifinal-Fig13-ses-dnn-mh.png", width = 3000, height = 1500)
op <- par(mfcol = c(2, 3), mar = c(4, 4, 0.5, 0.5), cex = 2)

plot.ses(beta.ses$tdnn0[c(T, F)], mh[c(T, F)], col = "tomato", lab = expression("SES " * D[nn]), xlab = "Mean altitude, m")
legend("topright", legend = "Tree layer", bty = "n")
plot.ses(beta.ses$tdnn1[c(T, F)], mh[c(T, F)], col = "tomato", lab = expression("SES D'"[nn]), xlab = "Mean altitude, m")

plot.ses(beta.ses$sdnn0[c(T, F)], mh[c(T, F)], col = "skyblue", lab = expression("SES " * D[nn]), xlab = "Mean altitude, m")
legend("topright", legend = "Shrub layer", bty = "n")
plot.ses(beta.ses$sdnn1[c(T, F)], mh[c(T, F)], col = "skyblue", lab = expression("SES D'"[nn]), xlab = "Mean altitude, m")

plot.ses(beta.ses$hdnn0[c(T, F)], mh[c(T, F)], col = "forestgreen", lab = expression("SES " * D[nn]), xlab = "Mean altitude, m")
legend("topright", legend = "Herb layer", bty = "n")
plot.ses(beta.ses$hdnn1[c(T, F)], mh[c(T, F)], col = "forestgreen", lab = expression("SES D'"[nn]), xlab = "Mean altitude, m")

par(op)
dev.off()

