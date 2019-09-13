library(vegan)
library(picante)

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

library(ecomf)

tpd <- PqD(ta, wood.tree)
spd <- PqD(sa, wood.tree)
hpd <- PqD(ha, herb.tree)

png("figures/qfinal-Fig01a.png", width = 1200, height = 800)
op <- par(mar = c(4, 4, 0.5, 0.5), cex = 2)

plot(hgt, tpd[, 1], pch = 19, col = "tomato", xlab = "Altitude, m", ylab = expression(PD[0]))
fit1 <- lm(tpd[, 1] ~ hgt)
fit2 <- lm(tpd[, 1] ~ hgt + I(hgt^2))
sm1 <- summary(fit1)
sm2 <- summary(fit2)
p1 <- pf(sm1$fstatistic[1], df1 = sm1$fstatistic[2], df2 = sm1$fstatistic[3], lower.tail = F)
p2 <- pf(sm2$fstatistic[1], df1 = sm2$fstatistic[2], df2 = sm2$fstatistic[3], lower.tail = F)

usr <- par("usr")

n <- length(tpd[, 1])
text(usr[2], usr[4] - 0.25*(usr[4] - usr[3]) / 10, pos = 2, adj = c(1, 1),
     labels = "Tree layer")

text(usr[2], usr[4] - 0.75*(usr[4] - usr[3]) / 10, pos = 2, adj = c(1, 1),
     labels = bquote("Linear fit: p = " ~ .(round(p1, digits = 3)) ~
                       ", " ~ r^2 ~ " = " ~ .(round(sm1$r.squared, digits = 3)) ~
                       ", AIC = " ~ .(round(AIC(fit1) + (2*3*4/(n-3-1)), digits = 3))))

text(usr[2], usr[4] - 1.25*(usr[4] - usr[3]) / 10, pos = 2, adj = c(1, 1),
     labels = bquote("Quadratic fit: p = " ~ .(round(p2, digits = 3)) ~
                           ", " ~ r^2 ~ " = " ~ .(round(sm2$r.squared, digits = 3)) ~
                           ", AIC = " ~ .(round(AIC(fit2) + (2*4*5/(n-4-1)), digits = 3))))

x <- seq(min(hgt), max(hgt), len = 1000)
abc <- coef(fit2)
y <- abc[1] + abc[2] * x + abc[3] * x^2
lines(x, y)

dev.off()

png("figures/qfinal-Fig01b.png", width = 1200, height = 800)
op <- par(mar = c(4, 4, 0.5, 0.5), cex = 2)

plot.ses(tpd[, 3], hgt, col = "tomato", lab = expression(PD[2]), xlab = "Altitude, m")
legend("topright", legend = "Tree layer", bty = "n")

dev.off()


png("figures/qfinal-Fig01c.png", width = 1200, height = 800)
op <- par(mar = c(4, 4, 0.5, 0.5), cex = 2)

plot.ses(spd[, 1], hgt, col = "skyblue", lab = expression(PD[0]), xlab = "Altitude, m")
legend("topright", legend = "Shrub layer", bty = "n")

dev.off()

png("figures/qfinal-Fig01d.png", width = 1200, height = 800)
op <- par(mar = c(4, 4, 0.5, 0.5), cex = 2)

plot.ses(spd[, 3], hgt, col = "skyblue", lab = expression(PD[2]), xlab = "Altitude, m")
legend("topright", legend = "Shrub layer", bty = "n")

dev.off()


png("figures/qfinal-Fig01e.png", width = 1200, height = 800)
op <- par(mar = c(4, 4, 0.5, 0.5), cex = 2)

plot.ses(hpd[, 1], hgt, col = "forestgreen", lab = expression(PD[0]), xlab = "Altitude, m")
legend("topleft", legend = "Herb layer", bty = "n")
dev.off()


png("figures/qfinal-Fig01f.png", width = 1200, height = 800)
op <- par(mar = c(4, 4, 0.5, 0.5), cex = 2)

plot.ses(hpd[, 3], hgt, col = "forestgreen", lab = expression(PD[2]), xlab = "Altitude, m")
legend("topright", legend = "Herb layer", bty = "n")

dev.off()



corTPalpha <- rbind(c(cor(tpd[, 1], specnumber(ta)),
                    cor(spd[, 1], specnumber(sa)),
                    cor(hpd[, 1], specnumber(ha))),
                    c(cor(tpd[, 3], diversity(ta, "invsimpson")),
                      cor(spd[, 3], diversity(sa, "invsimpson")),
                      cor(hpd[, 3], diversity(ha, "invsimpson"))))

#---#

load("clean.data/phylo-ses-pqd.rda")


png("figures/qfinal-Fig02a.png", width = 1200, height = 800)
op <- par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot.ses(ta.ses$mpd, hgt, col = "tomato", lab = "SES MPD", xlab = "Altitude, m")
legend("topright", legend = "Tree layer", bty = "n")
dev.off()


png("figures/qfinal-Fig02b.png", width = 1200, height = 800)
op <- par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot.ses(ta.ses$mpd.a, hgt, col = "tomato", lab = expression("SES " * MPD[a]), xlab = "Altitude, m")
legend("topright", legend = "Tree layer", bty = "n")
dev.off()

png("figures/qfinal-Fig02c.png", width = 1200, height = 800)
op <- par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot.ses(sa.ses$mpd, hgt, col = "skyblue", lab = "SES MPD", xlab = "Altitude, m")
legend("topright", legend = "Shrub layer", bty = "n")
dev.off()

png("figures/qfinal-Fig02d.png", width = 1200, height = 800)
op <- par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot.ses(sa.ses$mpd.a, hgt, col = "skyblue", lab = expression("SES " * MPD[a]), xlab = "Altitude, m")
legend("topright", legend = "Shrub layer", bty = "n")
dev.off()

png("figures/qfinal-Fig02e.png", width = 1200, height = 800)
op <- par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot.ses(ha.ses$mpd, hgt, col = "forestgreen", lab = "SES MPD", xlab = "Altitude, m")
legend("top", legend = "Herb layer", bty = "n")
dev.off()

png("figures/qfinal-Fig02f.png", width = 1200, height = 800)
op <- par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot.ses(ha.ses$mpd.a, hgt, col = "forestgreen", lab = expression("SES " * MPD[a]), xlab = "Altitude, m")
legend("top", legend = "Herb layer", bty = "n")
dev.off()


#---#

dh <- dist(hgt)

ta.sor <- vegdist(ta, method = "bray", binary = T)
ta.mh <- vegdist(ta, method = "horn", binary = F)
sa.sor <- vegdist(sa, method = "bray", binary = T)
sa.mh <- vegdist(sa, method = "horn", binary = F)
ha.sor <- vegdist(ha, method = "bray", binary = T)
ha.mh <- vegdist(ha, method = "horn", binary = F)


ta.sor.95 <- diag(as.matrix(ta.sor)[-1, ])
ta.mh.95 <- diag(as.matrix(ta.mh)[-1, ])
sa.sor.95 <- diag(as.matrix(sa.sor)[-1, ])
sa.mh.95 <- diag(as.matrix(sa.mh)[-1, ])
ha.sor.95 <- diag(as.matrix(ha.sor)[-1, ])
ha.mh.95 <- diag(as.matrix(ha.mh)[-1, ])

mh <- (hgt[1:95] + hgt[2:96])/2

#---#

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


png("figures/qfinal-Fig03a.png", width = 1200, height = 800)
op <- par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot.ses(ta.dpw0.95[c(T, F)], mh[c(T, F)], col = "tomato", lab = expression(D[pw]), xlab = "Mean altitude, m")
legend("topright", legend = "Tree layer", bty = "n")
dev.off()

png("figures/qfinal-Fig03b.png", width = 1200, height = 800)
op <- par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot.ses(ta.dpw1.95[c(T, F)], mh[c(T, F)], col = "tomato", lab = expression("D'"[pw]), xlab = "Mean altitude, m")
legend("topright", legend = "Tree layer", bty = "n")
dev.off()

png("figures/qfinal-Fig03c.png", width = 1200, height = 800)
op <- par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot.ses(sa.dpw0.95[c(T, F)], mh[c(T, F)], col = "skyblue", lab = expression(D[pw]), xlab = "Mean altitude, m")
legend("topright", legend = "Shrub layer", bty = "n")
dev.off()

png("figures/qfinal-Fig03d.png", width = 1200, height = 800)
op <- par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot.ses(sa.dpw1.95[c(T, F)], mh[c(T, F)], col = "skyblue", lab = expression("D'"[pw]), xlab = "Mean altitude, m")
legend("topright", legend = "Shrub layer", bty = "n")
dev.off()

png("figures/qfinal-Fig03e.png", width = 1200, height = 800)
op <- par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot.ses(ha.dpw0.95[c(T, F)], mh[c(T, F)], col = "forestgreen", lab = expression(D[pw]), xlab = "Mean altitude, m")
legend("topright", legend = "Herb layer", bty = "n")
dev.off()

png("figures/qfinal-Fig03f.png", width = 1200, height = 800)
op <- par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot.ses(ha.dpw1.95[c(T, F)], mh[c(T, F)], col = "forestgreen", lab = expression("D'"[pw]), xlab = "Mean altitude, m")
legend("topright", legend = "Herb layer", bty = "n")
dev.off()



corTPbetaPW <- rbind(c(cor(ta.dpw0.95[c(T, F)], ta.sor.95[c(T, F)]),
                      cor(sa.dpw0.95[c(T, F)], sa.sor.95[c(T, F)]),
                      cor(ha.dpw0.95[c(T, F)], ha.sor.95[c(T, F)])),
                    c(cor(ta.dpw1.95[c(T, F)], ta.mh.95[c(T, F)]),
                      cor(sa.dpw1.95[c(T, F)], sa.mh.95[c(T, F)]),
                      cor(ha.dpw1.95[c(T, F)], ha.mh.95[c(T, F)])))


corTPbetaPW.p <- rbind(c(cor.test(ta.dpw0.95[c(T, F)], ta.sor.95[c(T, F)])$p.value,
                       cor.test(sa.dpw0.95[c(T, F)], sa.sor.95[c(T, F)])$p.value,
                       cor.test(ha.dpw0.95[c(T, F)], ha.sor.95[c(T, F)])$p.value),
                     c(cor.test(ta.dpw1.95[c(T, F)], ta.mh.95[c(T, F)])$p.value,
                       cor.test(sa.dpw1.95[c(T, F)], sa.mh.95[c(T, F)])$p.value,
                       cor.test(ha.dpw1.95[c(T, F)], ha.mh.95[c(T, F)])$p.value))



corTPbetaNN <- rbind(c(cor(ta.dnn0.95[c(T, F)], ta.sor.95[c(T, F)]),
                       cor(sa.dnn0.95[c(T, F)], sa.sor.95[c(T, F)]),
                       cor(ha.dnn0.95[c(T, F)], ha.sor.95[c(T, F)])),
                     c(cor(ta.dnn1.95[c(T, F)], ta.mh.95[c(T, F)]),
                       cor(sa.dnn1.95[c(T, F)], sa.mh.95[c(T, F)]),
                       cor(ha.dnn1.95[c(T, F)], ha.mh.95[c(T, F)])))

corTPbetaNN.p <- rbind(c(cor.test(ta.dnn0.95[c(T, F)], ta.sor.95[c(T, F)])$p.value,
                         cor.test(sa.dnn0.95[c(T, F)], sa.sor.95[c(T, F)])$p.value,
                         cor.test(ha.dnn0.95[c(T, F)], ha.sor.95[c(T, F)])$p.value),
                       c(cor.test(ta.dnn1.95[c(T, F)], ta.mh.95[c(T, F)])$p.value,
                         cor.test(sa.dnn1.95[c(T, F)], sa.mh.95[c(T, F)])$p.value,
                         cor.test(ha.dnn1.95[c(T, F)], ha.mh.95[c(T, F)])$p.value))

round(corTPalpha, dig = 3)
round(corTPbetaPW, dig = 3)
round(corTPbetaNN, dig = 3)

corTPbetaPW.p < 0.05
corTPbetaNN.p < 0.05

#---#

load("clean.data/pw-beta-tl.rda")

png("figures/qfinal-Fig04a.png", width = 1200, height = 800)
op <- par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot.ses(beta.ses$td0[c(T, F)], mh[c(T, F)], col = "tomato", lab = expression("SES " * D[pw]), xlab = "Mean altitude, m")
legend("topright", legend = "Tree layer", bty = "n")
dev.off()

png("figures/qfinal-Fig04b.png", width = 1200, height = 800)
op <- par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot.ses(beta.ses$td1[c(T, F)], mh[c(T, F)], col = "tomato", lab = expression("SES D'"[pw]), xlab = "Mean altitude, m")
legend("topright", legend = "Tree layer", bty = "n")
dev.off()

png("figures/qfinal-Fig04c.png", width = 1200, height = 800)
op <- par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot.ses(beta.ses$sd0[c(T, F)], mh[c(T, F)], col = "skyblue", lab = expression("SES " * D[pw]), xlab = "Mean altitude, m")
legend("topright", legend = "Shrub layer", bty = "n")
dev.off()

png("figures/qfinal-Fig04d.png", width = 1200, height = 800)
op <- par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot.ses(beta.ses$sd1[c(T, F)], mh[c(T, F)], col = "skyblue", lab = expression("SES D'"[pw]), xlab = "Mean altitude, m")
legend("topright", legend = "Shrub layer", bty = "n")
dev.off()

png("figures/qfinal-Fig04e.png", width = 1200, height = 800)
op <- par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot.ses(beta.ses$hd0[c(T, F)], mh[c(T, F)], col = "forestgreen", lab = expression("SES " * D[pw]), xlab = "Mean altitude, m")
legend("topright", legend = "Herb layer", bty = "n")
dev.off()

png("figures/qfinal-Fig04f.png", width = 1200, height = 800)
op <- par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot.ses(beta.ses$hd1[c(T, F)], mh[c(T, F)], col = "forestgreen", lab = expression("SES D'"[pw]), xlab = "Mean altitude, m")
legend("topright", legend = "Herb layer", bty = "n")
dev.off()
