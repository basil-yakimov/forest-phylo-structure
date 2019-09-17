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

#plot.ses(tpd[, 3], hgt, col = "tomato", lab = expression(PD[2]), xlab = "Altitude, m")
#legend("topright", legend = "Tree layer", bty = "n")

plot(hgt, tpd[, 3], pch = 19, col = "tomato", xlab = "Altitude, m", ylab = expression(PD[2]))
fit1 <- lm(tpd[, 3] ~ hgt)
fit2 <- lm(tpd[, 3] ~ hgt + I(hgt^2))
sm1 <- summary(fit1)
sm2 <- summary(fit2)
p1 <- pf(sm1$fstatistic[1], df1 = sm1$fstatistic[2], df2 = sm1$fstatistic[3], lower.tail = F)
p2 <- pf(sm2$fstatistic[1], df1 = sm2$fstatistic[2], df2 = sm2$fstatistic[3], lower.tail = F)

usr <- par("usr")

n <- length(tpd[, 3])
text(usr[2], usr[4] - 0.1 - 0.25*(usr[4] - usr[3]) / 10, pos = 2, adj = c(1, 1),
     labels = "Tree layer")

text(usr[2], usr[4] - 0.1 - 0.75*(usr[4] - usr[3]) / 10, pos = 2, adj = c(1, 1),
     labels = bquote("Linear fit: p < 0.001" ~
                       ", " ~ r^2 ~ " = " ~ .(round(sm1$r.squared, digits = 3)) ~
                       ", AIC = " ~ .(round(AIC(fit1) + (2*3*4/(n-3-1)), digits = 3))))

text(usr[2], usr[4] - 0.1 - 1.25*(usr[4] - usr[3]) / 10, pos = 2, adj = c(1, 1),
     labels = bquote("Quadratic fit: p < 0.001" ~ 
                       ", " ~ r^2 ~ " = " ~ .(round(sm2$r.squared, digits = 3)) ~
                       ", AIC = " ~ .(round(AIC(fit2) + (2*4*5/(n-4-1)), digits = 3))))

x <- seq(min(hgt), max(hgt), len = 1000)
abc <- coef(fit2)
y <- abc[1] + abc[2] * x + abc[3] * x^2
lines(x, y)

dev.off()



png("figures/qfinal-Fig01c.png", width = 1200, height = 800)
op <- par(mar = c(4, 4, 0.5, 0.5), cex = 2)

#plot.ses(spd[, 1], hgt, col = "skyblue", lab = expression(PD[0]), xlab = "Altitude, m")
#legend("topright", legend = "Shrub layer", bty = "n")

plot(hgt, spd[, 1], pch = 19, col = "skyblue", xlab = "Altitude, m", ylab = expression(PD[0]))
fit1 <- lm(spd[, 1] ~ hgt)
fit2 <- lm(spd[, 1] ~ hgt + I(hgt^2))
sm1 <- summary(fit1)
sm2 <- summary(fit2)
p1 <- pf(sm1$fstatistic[1], df1 = sm1$fstatistic[2], df2 = sm1$fstatistic[3], lower.tail = F)
p2 <- pf(sm2$fstatistic[1], df1 = sm2$fstatistic[2], df2 = sm2$fstatistic[3], lower.tail = F)

usr <- par("usr")

n <- length(spd[, 1])
text(usr[2], usr[4] - 0.1 - 0.25*(usr[4] - usr[3]) / 10, pos = 2, adj = c(1, 1),
     labels = "Shrub layer")

text(usr[2], usr[4] - 0.1 - 0.75*(usr[4] - usr[3]) / 10, pos = 2, adj = c(1, 1),
     labels = bquote("Linear fit: p < 0.001" ~
                       ", " ~ r^2 ~ " = " ~ .(round(sm1$r.squared, digits = 3)) ~
                       ", AIC = " ~ .(round(AIC(fit1) + (2*3*4/(n-3-1)), digits = 3))))

text(usr[2], usr[4] - 0.1 - 1.25*(usr[4] - usr[3]) / 10, pos = 2, adj = c(1, 1),
     labels = bquote("Quadratic fit: p < 0.001" ~
                       ", " ~ r^2 ~ " = " ~ .(round(sm2$r.squared, digits = 3)) ~
                       ", AIC = " ~ .(round(AIC(fit2) + (2*4*5/(n-4-1)), digits = 3))))

abline(fit1)

dev.off()




png("figures/qfinal-Fig01d.png", width = 1200, height = 800)
op <- par(mar = c(4, 4, 0.5, 0.5), cex = 2)

#plot.ses(spd[, 3], hgt, col = "skyblue", lab = expression(PD[2]), xlab = "Altitude, m")
#legend("topright", legend = "Shrub layer", bty = "n")

plot(hgt, spd[, 3], pch = 19, col = "skyblue", xlab = "Altitude, m", ylab = expression(PD[2]))
fit1 <- lm(spd[, 3] ~ hgt)
fit2 <- lm(spd[, 3] ~ hgt + I(hgt^2))
sm1 <- summary(fit1)
sm2 <- summary(fit2)
p1 <- pf(sm1$fstatistic[1], df1 = sm1$fstatistic[2], df2 = sm1$fstatistic[3], lower.tail = F)
p2 <- pf(sm2$fstatistic[1], df1 = sm2$fstatistic[2], df2 = sm2$fstatistic[3], lower.tail = F)

usr <- par("usr")

n <- length(spd[, 3])
text(usr[2], usr[4] - 0.1 - 0.25*(usr[4] - usr[3]) / 10, pos = 2, adj = c(1, 1),
     labels = "Shrub layer")

text(usr[2], usr[4] - 0.1 - 0.75*(usr[4] - usr[3]) / 10, pos = 2, adj = c(1, 1),
     labels = bquote("Linear fit: p < 0.001" ~
                       ", " ~ r^2 ~ " = " ~ .(round(sm1$r.squared, digits = 3)) ~
                       ", AIC = " ~ .(round(AIC(fit1) + (2*3*4/(n-3-1)), digits = 3))))

text(usr[2], usr[4] - 0.1 - 1.25*(usr[4] - usr[3]) / 10, pos = 2, adj = c(1, 1),
     labels = bquote("Quadratic fit: p < 0.001" ~
                       ", " ~ r^2 ~ " = " ~ .(round(sm2$r.squared, digits = 3)) ~
                       ", AIC = " ~ .(round(AIC(fit2) + (2*4*5/(n-4-1)), digits = 3))))

abline(fit1)

dev.off()



png("figures/qfinal-Fig01e.png", width = 1200, height = 800)
op <- par(mar = c(4, 4, 0.5, 0.5), cex = 2)

#plot.ses(hpd[, 1], hgt, col = "forestgreen", lab = expression(PD[0]), xlab = "Altitude, m")
#legend("topleft", legend = "Herb layer", bty = "n")

plot(hgt, hpd[, 1], pch = 19, col = "forestgreen", xlab = "Altitude, m", ylab = expression(PD[0]))
fit1 <- lm(hpd[, 1] ~ hgt)
fit2 <- lm(hpd[, 1] ~ hgt + I(hgt^2))
sm1 <- summary(fit1)
sm2 <- summary(fit2)
p1 <- pf(sm1$fstatistic[1], df1 = sm1$fstatistic[2], df2 = sm1$fstatistic[3], lower.tail = F)
p2 <- pf(sm2$fstatistic[1], df1 = sm2$fstatistic[2], df2 = sm2$fstatistic[3], lower.tail = F)

usr <- par("usr")

n <- length(hpd[, 1])
text(usr[1], usr[4] - 0.1 - 0.25*(usr[4] - usr[3]) / 10, pos = 4, adj = c(1, 1),
     labels = "Herb layer")

text(usr[1], usr[4] - 0.1 - 0.75*(usr[4] - usr[3]) / 10, pos = 4, adj = c(1, 1),
     labels = bquote("Linear fit: p < 0.001" ~
                       ", " ~ r^2 ~ " = " ~ .(round(sm1$r.squared, digits = 3)) ~
                       ", AIC = " ~ .(round(AIC(fit1) + (2*3*4/(n-3-1)), digits = 3))))

text(usr[1], usr[4] - 0.1 - 1.25*(usr[4] - usr[3]) / 10, pos = 4, adj = c(1, 1),
     labels = bquote("Quadratic fit: p < 0.001" ~
                       ", " ~ r^2 ~ " = " ~ .(round(sm2$r.squared, digits = 3)) ~
                       ", AIC = " ~ .(round(AIC(fit2) + (2*4*5/(n-4-1)), digits = 3))))

x <- seq(min(hgt), max(hgt), len = 1000)
abc <- coef(fit2)
y <- abc[1] + abc[2] * x + abc[3] * x^2
lines(x, y)

dev.off()



png("figures/qfinal-Fig01f.png", width = 1200, height = 800)
op <- par(mar = c(4, 4, 0.5, 0.5), cex = 2)

#plot.ses(hpd[, 3], hgt, col = "forestgreen", lab = expression(PD[2]), xlab = "Altitude, m")
#legend("topright", legend = "Herb layer", bty = "n")

plot(hgt, hpd[, 3], pch = 19, col = "forestgreen", xlab = "Altitude, m", ylab = expression(PD[2]), ylim = c(1.334367, 4.4))
fit1 <- lm(hpd[, 3] ~ hgt)
fit2 <- lm(hpd[, 3] ~ hgt + I(hgt^2))
sm1 <- summary(fit1)
sm2 <- summary(fit2)
p1 <- pf(sm1$fstatistic[1], df1 = sm1$fstatistic[2], df2 = sm1$fstatistic[3], lower.tail = F)
p2 <- pf(sm2$fstatistic[1], df1 = sm2$fstatistic[2], df2 = sm2$fstatistic[3], lower.tail = F)

usr <- par("usr")

n <- length(hpd[, 3])
text(usr[2], usr[4] - 0.25*(usr[4] - usr[3]) / 10, pos = 2, adj = c(1, 1),
     labels = "Herb layer")

text(usr[2], usr[4] - 0.75*(usr[4] - usr[3]) / 10, pos = 2, adj = c(1, 1),
     labels = bquote("Linear fit: p < 0.001" ~
                       ", " ~ r^2 ~ " = " ~ .(round(sm1$r.squared, digits = 3)) ~
                       ", AIC = " ~ .(round(AIC(fit1) + (2*3*4/(n-3-1)), digits = 3))))

text(usr[2], usr[4] - 1.25*(usr[4] - usr[3]) / 10, pos = 2, adj = c(1, 1),
     labels = bquote("Quadratic fit: p < 0.001" ~
                       ", " ~ r^2 ~ " = " ~ .(round(sm2$r.squared, digits = 3)) ~
                       ", AIC = " ~ .(round(AIC(fit2) + (2*4*5/(n-4-1)), digits = 3))))

x <- seq(min(hgt), max(hgt), len = 1000)
abc <- coef(fit2)
y <- abc[1] + abc[2] * x + abc[3] * x^2
lines(x, y)

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

#plot.ses(ta.ses$mpd, hgt, col = "tomato", lab = "SES MPD", xlab = "Altitude, m")
#legend("topright", legend = "Tree layer", bty = "n")

plot(hgt, ta.ses$mpd, pch = 19, col = "tomato", xlab = "Altitude, m", ylab = "SES MPD")
fit1 <- lm(ta.ses$mpd ~ hgt)
fit2 <- lm(ta.ses$mpd ~ hgt + I(hgt^2))
sm1 <- summary(fit1)
sm2 <- summary(fit2)
p1 <- pf(sm1$fstatistic[1], df1 = sm1$fstatistic[2], df2 = sm1$fstatistic[3], lower.tail = F)
p2 <- pf(sm2$fstatistic[1], df1 = sm2$fstatistic[2], df2 = sm2$fstatistic[3], lower.tail = F)

usr <- par("usr")

n <- length(ta.ses$mpd)
text(usr[1], usr[3] + 0.1 - 1.25*(usr[3] - usr[4]) / 10, pos = 4, adj = c(1, 1),
     labels = "Tree layer")

text(usr[1], usr[3] + 0.1 - 0.75*(usr[3] - usr[4]) / 10, pos = 4, adj = c(1, 1),
     labels = bquote("Linear fit: p < 0.001" ~
                       ", " ~ r^2 ~ " = " ~ .(round(sm1$r.squared, digits = 3)) ~
                       ", AIC = " ~ .(round(AIC(fit1) + (2*3*4/(n-3-1)), digits = 3))))

text(usr[1], usr[3] + 0.1 - 0.25*(usr[3] - usr[4]) / 10, pos = 4, adj = c(1, 1),
     labels = bquote("Quadratic fit:  p = " ~ .(round(p2, digits = 3)) ~
                       ", " ~ r^2 ~ " = " ~ .(round(sm2$r.squared, digits = 3)) ~
                       ", AIC = " ~ .(round(AIC(fit2) + (2*4*5/(n-4-1)), digits = 3))))

abline(fit1)

dev.off()


png("figures/qfinal-Fig02b.png", width = 1200, height = 800)
op <- par(mar = c(4, 4, 0.5, 0.5), cex = 2)
#plot.ses(ta.ses$mpd.a, hgt, col = "tomato", lab = expression("SES " * MPD[a]), xlab = "Altitude, m")
#legend("topright", legend = "Tree layer", bty = "n")

plot(hgt, ta.ses$mpd.a, pch = 19, col = "tomato", xlab = "Altitude, m", ylab = expression("SES " * MPD[a]))
fit1 <- lm(ta.ses$mpd.a ~ hgt)
fit2 <- lm(ta.ses$mpd.a ~ hgt + I(hgt^2))
sm1 <- summary(fit1)
sm2 <- summary(fit2)
p1 <- pf(sm1$fstatistic[1], df1 = sm1$fstatistic[2], df2 = sm1$fstatistic[3], lower.tail = F)
p2 <- pf(sm2$fstatistic[1], df1 = sm2$fstatistic[2], df2 = sm2$fstatistic[3], lower.tail = F)

usr <- par("usr")

n <- length(ta.ses$mpd.a)
text(usr[2], usr[3] + 0.1 - 1.25*(usr[3] - usr[4]) / 10, pos = 2, adj = c(1, 1),
     labels = "Tree layer")

text(usr[2], usr[3] + 0.1 - 0.75*(usr[3] - usr[4]) / 10, pos = 2, adj = c(1, 1),
     labels = bquote("Linear fit: p < 0.001" ~
                       ", " ~ r^2 ~ " = " ~ .(round(sm1$r.squared, digits = 3)) ~
                       ", AIC = " ~ .(round(AIC(fit1) + (2*3*4/(n-3-1)), digits = 3))))

text(usr[2], usr[3] + 0.1 - 0.25*(usr[3] - usr[4]) / 10, pos = 2, adj = c(1, 1),
     labels = bquote("Quadratic fit:  p < 0.001" ~
                       ", " ~ r^2 ~ " = " ~ .(round(sm2$r.squared, digits = 3)) ~
                       ", AIC = " ~ .(round(AIC(fit2) + (2*4*5/(n-4-1)), digits = 3))))

abline(fit1)

dev.off()

png("figures/qfinal-Fig02c.png", width = 1200, height = 800)
op <- par(mar = c(4, 4, 0.5, 0.5), cex = 2)
#plot.ses(sa.ses$mpd, hgt, col = "skyblue", lab = "SES MPD", xlab = "Altitude, m")
#legend("topright", legend = "Shrub layer", bty = "n")

plot(hgt, sa.ses$mpd, pch = 19, col = "skyblue", xlab = "Altitude, m", ylab = "SES MPD", ylim = c(-3.925523, 2.15))
fit1 <- lm(sa.ses$mpd ~ hgt)
fit2 <- lm(sa.ses$mpd ~ hgt + I(hgt^2))
sm1 <- summary(fit1)
sm2 <- summary(fit2)
p1 <- pf(sm1$fstatistic[1], df1 = sm1$fstatistic[2], df2 = sm1$fstatistic[3], lower.tail = F)
p2 <- pf(sm2$fstatistic[1], df1 = sm2$fstatistic[2], df2 = sm2$fstatistic[3], lower.tail = F)

usr <- par("usr")

n <- length(sa.ses$mpd)
text(usr[2], usr[4] - 0.25*(usr[4] - usr[3]) / 10, pos = 2, adj = c(1, 1),
     labels = "Shrub layer")

text(usr[2], usr[4] - 0.75*(usr[4] - usr[3]) / 10, pos = 2, adj = c(1, 1),
     labels = bquote("Linear fit: p < 0.001" ~
                       ", " ~ r^2 ~ " = " ~ .(round(sm1$r.squared, digits = 3)) ~
                       ", AIC = " ~ .(round(AIC(fit1) + (2*3*4/(n-3-1)), digits = 3))))

text(usr[2], usr[4] - 1.25*(usr[4] - usr[3]) / 10, pos = 2, adj = c(1, 1),
     labels = bquote("Quadratic fit: p < 0.001" ~
                       ", " ~ r^2 ~ " = " ~ .(round(sm2$r.squared, digits = 3)) ~
                       ", AIC = " ~ .(round(AIC(fit2) + (2*4*5/(n-4-1)), digits = 3))))

abline(fit1)

dev.off()

png("figures/qfinal-Fig02d.png", width = 1200, height = 800)
op <- par(mar = c(4, 4, 0.5, 0.5), cex = 2)
#plot.ses(sa.ses$mpd.a, hgt, col = "skyblue", lab = expression("SES " * MPD[a]), xlab = "Altitude, m")
#legend("topright", legend = "Shrub layer", bty = "n")

plot(hgt, sa.ses$mpd.a, pch = 19, col = "skyblue", xlab = "Altitude, m", ylab = expression("SES " * MPD[a]), ylim = c(-1.587927, 2))
fit1 <- lm(sa.ses$mpd.a ~ hgt)
fit2 <- lm(sa.ses$mpd.a ~ hgt + I(hgt^2))
sm1 <- summary(fit1)
sm2 <- summary(fit2)
p1 <- pf(sm1$fstatistic[1], df1 = sm1$fstatistic[2], df2 = sm1$fstatistic[3], lower.tail = F)
p2 <- pf(sm2$fstatistic[1], df1 = sm2$fstatistic[2], df2 = sm2$fstatistic[3], lower.tail = F)

usr <- par("usr")

n <- length(sa.ses$mpd.a)
text(usr[2], usr[4] - 0.25*(usr[4] - usr[3]) / 10, pos = 2, adj = c(1, 1),
     labels = "Shrub layer")

text(usr[2], usr[4] - 0.75*(usr[4] - usr[3]) / 10, pos = 2, adj = c(1, 1),
     labels = bquote("Linear fit: p = " ~ .(round(p1, digits = 3)) ~
                       ", " ~ r^2 ~ " = " ~ .(round(sm1$r.squared, digits = 3)) ~
                       ", AIC = " ~ .(round(AIC(fit1) + (2*3*4/(n-3-1)), digits = 3))))

text(usr[2], usr[4] - 1.25*(usr[4] - usr[3]) / 10, pos = 2, adj = c(1, 1),
     labels = bquote("Quadratic fit:  p = " ~ .(round(p2, digits = 3)) ~
                       ", " ~ r^2 ~ " = " ~ .(round(sm2$r.squared, digits = 3)) ~
                       ", AIC = " ~ .(round(AIC(fit2) + (2*4*5/(n-4-1)), digits = 3))))

x <- seq(min(hgt), max(hgt), len = 1000)
abc <- coef(fit2)
y <- abc[1] + abc[2] * x + abc[3] * x^2
lines(x, y)

dev.off()

png("figures/qfinal-Fig02e.png", width = 1200, height = 800)
op <- par(mar = c(4, 4, 0.5, 0.5), cex = 2)
#plot.ses(ha.ses$mpd, hgt, col = "forestgreen", lab = "SES MPD", xlab = "Altitude, m")
#legend("top", legend = "Herb layer", bty = "n")

plot(hgt, ha.ses$mpd, pch = 19, col = "forestgreen", xlab = "Altitude, m", ylab = "SES MPD", ylim = c(-1.587927, 2.2))
fit1 <- lm(ha.ses$mpd ~ hgt)
fit2 <- lm(ha.ses$mpd ~ hgt + I(hgt^2))
sm1 <- summary(fit1)
sm2 <- summary(fit2)
p1 <- pf(sm1$fstatistic[1], df1 = sm1$fstatistic[2], df2 = sm1$fstatistic[3], lower.tail = F)
p2 <- pf(sm2$fstatistic[1], df1 = sm2$fstatistic[2], df2 = sm2$fstatistic[3], lower.tail = F)

usr <- par("usr")

n <- length(ha.ses$mpd)
text(usr[1], usr[4] - 0.25*(usr[4] - usr[3]) / 10, pos = 4, adj = c(1, 1),
     labels = "Herb layer")

text(usr[1], usr[4] - 0.75*(usr[4] - usr[3]) / 10, pos = 4, adj = c(1, 1),
     labels = bquote("Linear fit: p = " ~ .(round(p1, digits = 3)) ~
                       ", " ~ r^2 ~ " = " ~ .(round(sm1$r.squared, digits = 3)) ~
                       ", AIC = " ~ .(round(AIC(fit1) + (2*3*4/(n-3-1)), digits = 3))))

text(usr[1], usr[4] - 1.25*(usr[4] - usr[3]) / 10, pos = 4, adj = c(1, 1),
     labels = bquote("Quadratic fit:  p = " ~ .(round(p2, digits = 3)) ~
                       ", " ~ r^2 ~ " = " ~ .(round(sm2$r.squared, digits = 3)) ~
                       ", AIC = " ~ .(round(AIC(fit2) + (2*4*5/(n-4-1)), digits = 3))))

dev.off()

png("figures/qfinal-Fig02f.png", width = 1200, height = 800)
op <- par(mar = c(4, 4, 0.5, 0.5), cex = 2)
#plot.ses(ha.ses$mpd.a, hgt, col = "forestgreen", lab = expression("SES " * MPD[a]), xlab = "Altitude, m")
#legend("top", legend = "Herb layer", bty = "n")

plot(hgt, ha.ses$mpd.a, pch = 19, col = "forestgreen", xlab = "Altitude, m", ylab = expression("SES " * MPD[a]), ylim = c(-2.270101, 2))
fit1 <- lm(ha.ses$mpd.a ~ hgt)
fit2 <- lm(ha.ses$mpd.a ~ hgt + I(hgt^2))
sm1 <- summary(fit1)
sm2 <- summary(fit2)
p1 <- pf(sm1$fstatistic[1], df1 = sm1$fstatistic[2], df2 = sm1$fstatistic[3], lower.tail = F)
p2 <- pf(sm2$fstatistic[1], df1 = sm2$fstatistic[2], df2 = sm2$fstatistic[3], lower.tail = F)

usr <- par("usr")

n <- length(ha.ses$mpd.a)
text(usr[2], usr[4] - 0.25*(usr[4] - usr[3]) / 10, pos = 2, adj = c(1, 1),
     labels = "Herb layer")

text(usr[2], usr[4] - 0.75*(usr[4] - usr[3]) / 10, pos = 2, adj = c(1, 1),
     labels = bquote("Linear fit: p = " ~ .(round(p1, digits = 3)) ~
                       ", " ~ r^2 ~ " = " ~ .(round(sm1$r.squared, digits = 3)) ~
                       ", AIC = " ~ .(round(AIC(fit1) + (2*3*4/(n-3-1)), digits = 3))))

text(usr[2], usr[4] - 1.25*(usr[4] - usr[3]) / 10, pos = 2, adj = c(1, 1),
     labels = bquote("Quadratic fit:  p = " ~ .(round(p2, digits = 3)) ~
                       ", " ~ r^2 ~ " = " ~ .(round(sm2$r.squared, digits = 3)) ~
                       ", AIC = " ~ .(round(AIC(fit2) + (2*4*5/(n-4-1)), digits = 3))))
x <- seq(min(hgt), max(hgt), len = 1000)
abc <- coef(fit2)
y <- abc[1] + abc[2] * x + abc[3] * x^2
lines(x, y)

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
#plot.ses(ta.dpw0.95[c(T, F)], mh[c(T, F)], col = "tomato", lab = expression(D[pw]), xlab = "Mean altitude, m")
#legend("topright", legend = "Tree layer", bty = "n")

plot(mh[c(T, F)], ta.dpw0.95[c(T, F)], pch = 19, col = "tomato", xlab = "Mean altitude, m", ylab = expression(D[pw]))
fit1 <- lm(ta.dpw0.95[c(T, F)] ~ mh[c(T, F)])
fit2 <- lm(ta.dpw0.95[c(T, F)] ~ mh[c(T, F)] + I(mh[c(T, F)]^2))
sm1 <- summary(fit1)
sm2 <- summary(fit2)
p1 <- pf(sm1$fstatistic[1], df1 = sm1$fstatistic[2], df2 = sm1$fstatistic[3], lower.tail = F)
p2 <- pf(sm2$fstatistic[1], df1 = sm2$fstatistic[2], df2 = sm2$fstatistic[3], lower.tail = F)

usr <- par("usr")

n <- length(ta.dpw0.95[c(T, F)])
text(usr[1], usr[3] + 0.1 - 1.25*(usr[3] - usr[4]) / 10, pos = 4, adj = c(1, 1),
     labels = "Tree layer")

text(usr[1], usr[3] + 0.1 - 0.75*(usr[3] - usr[4]) / 10, pos = 4, adj = c(1, 1),
     labels = bquote("Linear fit: p < 0.001" ~
                       ", " ~ r^2 ~ " = " ~ .(round(sm1$r.squared, digits = 3)) ~
                       ", AIC = " ~ .(round(AIC(fit1) + (2*3*4/(n-3-1)), digits = 3))))

text(usr[1], usr[3] + 0.1 - 0.25*(usr[3] - usr[4]) / 10, pos = 4, adj = c(1, 1),
     labels = bquote("Quadratic fit:  p < 0.001" ~
                       ", " ~ r^2 ~ " = " ~ .(round(sm2$r.squared, digits = 3)) ~
                       ", AIC = " ~ .(round(AIC(fit2) + (2*4*5/(n-4-1)), digits = 3))))
x <- seq(min(hgt), max(hgt), len = 1000)
abc <- coef(fit2)
y <- abc[1] + abc[2] * x + abc[3] * x^2
lines(x, y)

dev.off()

png("figures/qfinal-Fig03b.png", width = 1200, height = 800)
op <- par(mar = c(4, 4, 0.5, 0.5), cex = 2)
#plot.ses(ta.dpw1.95[c(T, F)], mh[c(T, F)], col = "tomato", lab = expression("D'"[pw]), xlab = "Mean altitude, m")
#legend("topright", legend = "Tree layer", bty = "n")

plot(mh[c(T, F)], ta.dpw1.95[c(T, F)], pch = 19, col = "tomato", xlab = "Mean altitude, m", ylab = expression("D'"[pw]))
fit1 <- lm(ta.dpw1.95[c(T, F)] ~ mh[c(T, F)])
fit2 <- lm(ta.dpw1.95[c(T, F)] ~ mh[c(T, F)] + I(mh[c(T, F)]^2))
sm1 <- summary(fit1)
sm2 <- summary(fit2)
p1 <- pf(sm1$fstatistic[1], df1 = sm1$fstatistic[2], df2 = sm1$fstatistic[3], lower.tail = F)
p2 <- pf(sm2$fstatistic[1], df1 = sm2$fstatistic[2], df2 = sm2$fstatistic[3], lower.tail = F)

usr <- par("usr")

n <- length(ta.dpw1.95[c(T, F)])
text(usr[2], usr[4] - 0.1 - 0.25*(usr[4] - usr[3]) / 10, pos = 2, adj = c(1, 1),
     labels = "Tree layer")

text(usr[2], usr[4] - 0.1 - 0.75*(usr[4] - usr[3]) / 10, pos = 2, adj = c(1, 1),
     labels = bquote("Linear fit: p < 0.001" ~
                       ", " ~ r^2 ~ " = " ~ .(round(sm1$r.squared, digits = 3)) ~
                       ", AIC = " ~ .(round(AIC(fit1) + (2*3*4/(n-3-1)), digits = 3))))

text(usr[2], usr[4] - 0.1 - 1.25*(usr[4] - usr[3]) / 10, pos = 2, adj = c(1, 1),
     labels = bquote("Quadratic fit: p < 0.001" ~ 
                       ", " ~ r^2 ~ " = " ~ .(round(sm2$r.squared, digits = 3)) ~
                       ", AIC = " ~ .(round(AIC(fit2) + (2*4*5/(n-4-1)), digits = 3))))

x <- seq(min(hgt), max(hgt), len = 1000)
abc <- coef(fit2)
y <- abc[1] + abc[2] * x + abc[3] * x^2
lines(x, y)

dev.off()

png("figures/qfinal-Fig03c.png", width = 1200, height = 800)
op <- par(mar = c(4, 4, 0.5, 0.5), cex = 2)
#plot.ses(sa.dpw0.95[c(T, F)], mh[c(T, F)], col = "skyblue", lab = expression(D[pw]), xlab = "Mean altitude, m")
#legend("topright", legend = "Shrub layer", bty = "n")

plot(mh[c(T, F)], sa.dpw0.95[c(T, F)], pch = 19, col = "skyblue", xlab = "Mean altitude, m", ylab = expression(D[pw]))
fit1 <- lm(sa.dpw0.95[c(T, F)] ~ mh[c(T, F)])
fit2 <- lm(sa.dpw0.95[c(T, F)] ~ mh[c(T, F)] + I(mh[c(T, F)]^2))
sm1 <- summary(fit1)
sm2 <- summary(fit2)
p1 <- pf(sm1$fstatistic[1], df1 = sm1$fstatistic[2], df2 = sm1$fstatistic[3], lower.tail = F)
p2 <- pf(sm2$fstatistic[1], df1 = sm2$fstatistic[2], df2 = sm2$fstatistic[3], lower.tail = F)

usr <- par("usr")

n <- length(sa.dpw0.95[c(T, F)])
text(usr[1], usr[3] + 0.1 - 1.25*(usr[3] - usr[4]) / 10, pos = 4, adj = c(1, 1),
     labels = "Shrub layer")

text(usr[1], usr[3] + 0.1 - 0.75*(usr[3] - usr[4]) / 10, pos = 4, adj = c(1, 1),
     labels = bquote("Linear fit: p < 0.001" ~
                       ", " ~ r^2 ~ " = " ~ .(round(sm1$r.squared, digits = 3)) ~
                       ", AIC = " ~ .(round(AIC(fit1) + (2*3*4/(n-3-1)), digits = 3))))

text(usr[1], usr[3] + 0.1 - 0.25*(usr[3] - usr[4]) / 10, pos = 4, adj = c(1, 1),
     labels = bquote("Quadratic fit:  p < 0.001" ~
                       ", " ~ r^2 ~ " = " ~ .(round(sm2$r.squared, digits = 3)) ~
                       ", AIC = " ~ .(round(AIC(fit2) + (2*4*5/(n-4-1)), digits = 3))))
x <- seq(min(hgt), max(hgt), len = 1000)
abc <- coef(fit2)
y <- abc[1] + abc[2] * x + abc[3] * x^2
lines(x, y)

dev.off()

png("figures/qfinal-Fig03d.png", width = 1200, height = 800)
op <- par(mar = c(4, 4, 0.5, 0.5), cex = 2)
#plot.ses(sa.dpw1.95[c(T, F)], mh[c(T, F)], col = "skyblue", lab = expression("D'"[pw]), xlab = "Mean altitude, m")
#legend("topright", legend = "Shrub layer", bty = "n")

plot(mh[c(T, F)], sa.dpw1.95[c(T, F)], pch = 19, col = "skyblue", xlab = "Mean altitude, m", ylab = expression("D'"[pw]))
fit1 <- lm(sa.dpw1.95[c(T, F)] ~ mh[c(T, F)])
fit2 <- lm(sa.dpw1.95[c(T, F)] ~ mh[c(T, F)] + I(mh[c(T, F)]^2))
sm1 <- summary(fit1)
sm2 <- summary(fit2)
p1 <- pf(sm1$fstatistic[1], df1 = sm1$fstatistic[2], df2 = sm1$fstatistic[3], lower.tail = F)
p2 <- pf(sm2$fstatistic[1], df1 = sm2$fstatistic[2], df2 = sm2$fstatistic[3], lower.tail = F)

usr <- par("usr")

n <- length(sa.dpw1.95[c(T, F)])
text(usr[1], usr[3] + 0.1 - 1.25*(usr[3] - usr[4]) / 10, pos = 4, adj = c(1, 1),
     labels = "Shrub layer")

text(usr[1], usr[3] + 0.1 - 0.75*(usr[3] - usr[4]) / 10, pos = 4, adj = c(1, 1),
     labels = bquote("Linear fit: p < 0.001" ~
                       ", " ~ r^2 ~ " = " ~ .(round(sm1$r.squared, digits = 3)) ~
                       ", AIC = " ~ .(round(AIC(fit1) + (2*3*4/(n-3-1)), digits = 3))))

text(usr[1], usr[3] + 0.1 - 0.25*(usr[3] - usr[4]) / 10, pos = 4, adj = c(1, 1),
     labels = bquote("Quadratic fit:  p < 0.001" ~
                       ", " ~ r^2 ~ " = " ~ .(round(sm2$r.squared, digits = 3)) ~
                       ", AIC = " ~ .(round(AIC(fit2) + (2*4*5/(n-4-1)), digits = 3))))
x <- seq(min(hgt), max(hgt), len = 1000)
abc <- coef(fit2)
y <- abc[1] + abc[2] * x + abc[3] * x^2
lines(x, y)

dev.off()

png("figures/qfinal-Fig03e.png", width = 1200, height = 800)
op <- par(mar = c(4, 4, 0.5, 0.5), cex = 2)
#plot.ses(ha.dpw0.95[c(T, F)], mh[c(T, F)], col = "forestgreen", lab = expression(D[pw]), xlab = "Mean altitude, m")
#legend("topright", legend = "Herb layer", bty = "n")

plot(mh[c(T, F)], ha.dpw0.95[c(T, F)], pch = 19, col = "forestgreen", xlab = "Mean altitude, m", ylab = expression(D[pw]), ylim = c(233.5425, 308))
fit1 <- lm(ha.dpw0.95[c(T, F)] ~ mh[c(T, F)])
fit2 <- lm(ha.dpw0.95[c(T, F)] ~ mh[c(T, F)] + I(mh[c(T, F)]^2))
sm1 <- summary(fit1)
sm2 <- summary(fit2)
p1 <- pf(sm1$fstatistic[1], df1 = sm1$fstatistic[2], df2 = sm1$fstatistic[3], lower.tail = F)
p2 <- pf(sm2$fstatistic[1], df1 = sm2$fstatistic[2], df2 = sm2$fstatistic[3], lower.tail = F)

usr <- par("usr")

n <- length(ha.dpw0.95[c(T, F)])
text(usr[1], usr[4] - 0.25*(usr[4] - usr[3]) / 10, pos = 4, adj = c(1, 1),
     labels = "Herb layer")

text(usr[1], usr[4] - 0.75*(usr[4] - usr[3]) / 10, pos = 4, adj = c(1, 1),
     labels = bquote("Linear fit: p = " ~ .(round(p1, digits = 3)) ~
                       ", " ~ r^2 ~ " = " ~ .(round(sm1$r.squared, digits = 3)) ~
                       ", AIC = " ~ .(round(AIC(fit1) + (2*3*4/(n-3-1)), digits = 3))))

text(usr[1], usr[4] - 1.25*(usr[4] - usr[3]) / 10, pos = 4, adj = c(1, 1),
     labels = bquote("Quadratic fit:  p = " ~ .(round(p2, digits = 3)) ~
                       ", " ~ r^2 ~ " = " ~ .(round(sm2$r.squared, digits = 3)) ~
                       ", AIC = " ~ .(round(AIC(fit2) + (2*4*5/(n-4-1)), digits = 3))))

dev.off()

png("figures/qfinal-Fig03f.png", width = 1200, height = 800)
op <- par(mar = c(4, 4, 0.5, 0.5), cex = 2)
#plot.ses(ha.dpw1.95[c(T, F)], mh[c(T, F)], col = "forestgreen", lab = expression("D'"[pw]), xlab = "Mean altitude, m")
#legend("topright", legend = "Herb layer", bty = "n")

plot(mh[c(T, F)], ha.dpw1.95[c(T, F)], pch = 19, col = "forestgreen", xlab = "Mean altitude, m", ylab = expression("D'"[pw]), xlim = c(1060, 1750))
fit1 <- lm(ha.dpw1.95[c(T, F)] ~ mh[c(T, F)])
fit2 <- lm(ha.dpw1.95[c(T, F)] ~ mh[c(T, F)] + I(mh[c(T, F)]^2))
sm1 <- summary(fit1)
sm2 <- summary(fit2)
p1 <- pf(sm1$fstatistic[1], df1 = sm1$fstatistic[2], df2 = sm1$fstatistic[3], lower.tail = F)
p2 <- pf(sm2$fstatistic[1], df1 = sm2$fstatistic[2], df2 = sm2$fstatistic[3], lower.tail = F)

usr <- par("usr")

n <- length(ha.dpw1.95[c(T, F)])
text(usr[1], usr[3] + 0.1 - 1.25*(usr[3] - usr[4]) / 10, pos = 4, adj = c(1, 1),
     labels = "Herb layer")

text(usr[1], usr[3] + 0.1 - 0.75*(usr[3] - usr[4]) / 10, pos = 4, adj = c(1, 1),
     labels = bquote("Linear fit: p < 0.001" ~
                       ", " ~ r^2 ~ " = " ~ .(round(sm1$r.squared, digits = 3)) ~
                       ", AIC = " ~ .(round(AIC(fit1) + (2*3*4/(n-3-1)), digits = 3))))

text(usr[1], usr[3] + 0.1 - 0.25*(usr[3] - usr[4]) / 10, pos = 4, adj = c(1, 1),
     labels = bquote("Quadratic fit:  p < 0.001" ~
                       ", " ~ r^2 ~ " = " ~ .(round(sm2$r.squared, digits = 3)) ~
                       ", AIC = " ~ .(round(AIC(fit2) + (2*4*5/(n-4-1)), digits = 3))))
x <- seq(min(hgt), max(hgt), len = 1000)
abc <- coef(fit2)
y <- abc[1] + abc[2] * x + abc[3] * x^2
lines(x, y)

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
#plot.ses(beta.ses$td0[c(T, F)], mh[c(T, F)], col = "tomato", lab = expression("SES " * D[pw]), xlab = "Mean altitude, m")
#legend("topright", legend = "Tree layer", bty = "n")

plot(mh[c(T, F)], beta.ses$td0[c(T, F)], pch = 19, col = "tomato", xlab = "Mean altitude, m", ylab = expression("SES " * D[pw]), xlim = c(1060, 1750))
fit1 <- lm(beta.ses$td0[c(T, F)] ~ mh[c(T, F)])
fit2 <- lm(beta.ses$td0[c(T, F)] ~ mh[c(T, F)] + I(mh[c(T, F)]^2))
sm1 <- summary(fit1)
sm2 <- summary(fit2)
p1 <- pf(sm1$fstatistic[1], df1 = sm1$fstatistic[2], df2 = sm1$fstatistic[3], lower.tail = F)
p2 <- pf(sm2$fstatistic[1], df1 = sm2$fstatistic[2], df2 = sm2$fstatistic[3], lower.tail = F)

usr <- par("usr")

n <- length(beta.ses$td0[c(T, F)])
text(usr[1], usr[3] + 0.1 - 1.25*(usr[3] - usr[4]) / 10, pos = 4, adj = c(1, 1),
     labels = "Tree layer")

text(usr[1], usr[3] + 0.1 - 0.75*(usr[3] - usr[4]) / 10, pos = 4, adj = c(1, 1),
     labels = bquote("Linear fit: p < 0.001" ~
                       ", " ~ r^2 ~ " = " ~ .(round(sm1$r.squared, digits = 3)) ~
                       ", AIC = " ~ .(round(AIC(fit1) + (2*3*4/(n-3-1)), digits = 3))))

text(usr[1], usr[3] + 0.1 - 0.25*(usr[3] - usr[4]) / 10, pos = 4, adj = c(1, 1),
     labels = bquote("Quadratic fit:  p < 0.001" ~
                       ", " ~ r^2 ~ " = " ~ .(round(sm2$r.squared, digits = 3)) ~
                       ", AIC = " ~ .(round(AIC(fit2) + (2*4*5/(n-4-1)), digits = 3))))

abline(fit1)

dev.off()

png("figures/qfinal-Fig04b.png", width = 1200, height = 800)
op <- par(mar = c(4, 4, 0.5, 0.5), cex = 2)
#plot.ses(beta.ses$td1[c(T, F)], mh[c(T, F)], col = "tomato", lab = expression("SES D'"[pw]), xlab = "Mean altitude, m")
#legend("topright", legend = "Tree layer", bty = "n")

plot(mh[c(T, F)], beta.ses$td1[c(T, F)], pch = 19, col = "tomato", xlab = "Mean altitude, m", ylab = expression("SES " * D[pw]))
fit1 <- lm(beta.ses$td1[c(T, F)] ~ mh[c(T, F)])
fit2 <- lm(beta.ses$td1[c(T, F)] ~ mh[c(T, F)] + I(mh[c(T, F)]^2))
sm1 <- summary(fit1)
sm2 <- summary(fit2)
p1 <- pf(sm1$fstatistic[1], df1 = sm1$fstatistic[2], df2 = sm1$fstatistic[3], lower.tail = F)
p2 <- pf(sm2$fstatistic[1], df1 = sm2$fstatistic[2], df2 = sm2$fstatistic[3], lower.tail = F)

usr <- par("usr")

n <- length(beta.ses$td1[c(T, F)])
text(usr[1], usr[3] + 0.1 - 1.25*(usr[3] - usr[4]) / 10, pos = 4, adj = c(1, 1),
     labels = "Tree layer")

text(usr[1], usr[3] + 0.1 - 0.75*(usr[3] - usr[4]) / 10, pos = 4, adj = c(1, 1),
     labels = bquote("Linear fit: p < 0.001" ~
                       ", " ~ r^2 ~ " = " ~ .(round(sm1$r.squared, digits = 3)) ~
                       ", AIC = " ~ .(round(AIC(fit1) + (2*3*4/(n-3-1)), digits = 3))))

text(usr[1], usr[3] + 0.1 - 0.25*(usr[3] - usr[4]) / 10, pos = 4, adj = c(1, 1),
     labels = bquote("Quadratic fit:  p < 0.001" ~
                       ", " ~ r^2 ~ " = " ~ .(round(sm2$r.squared, digits = 3)) ~
                       ", AIC = " ~ .(round(AIC(fit2) + (2*4*5/(n-4-1)), digits = 3))))

abline(fit1)

dev.off()

png("figures/qfinal-Fig04c.png", width = 1200, height = 800)
op <- par(mar = c(4, 4, 0.5, 0.5), cex = 2)
#plot.ses(beta.ses$sd0[c(T, F)], mh[c(T, F)], col = "skyblue", lab = expression("SES " * D[pw]), xlab = "Mean altitude, m")
#legend("topright", legend = "Shrub layer", bty = "n")

plot(mh[c(T, F)], beta.ses$sd0[c(T, F)], pch = 19, col = "skyblue", xlab = "Mean altitude, m", ylab = expression("SES " * D[pw]))
fit1 <- lm(beta.ses$sd0[c(T, F)] ~ mh[c(T, F)])
fit2 <- lm(beta.ses$sd0[c(T, F)] ~ mh[c(T, F)] + I(mh[c(T, F)]^2))
sm1 <- summary(fit1)
sm2 <- summary(fit2)
p1 <- pf(sm1$fstatistic[1], df1 = sm1$fstatistic[2], df2 = sm1$fstatistic[3], lower.tail = F)
p2 <- pf(sm2$fstatistic[1], df1 = sm2$fstatistic[2], df2 = sm2$fstatistic[3], lower.tail = F)

usr <- par("usr")

n <- length(beta.ses$sd0[c(T, F)])
text(usr[1], usr[3] + 0.1 - 1.25*(usr[3] - usr[4]) / 10, pos = 4, adj = c(1, 1),
     labels = "Shrub layer")

text(usr[1], usr[3] + 0.1 - 0.75*(usr[3] - usr[4]) / 10, pos = 4, adj = c(1, 1),
     labels = bquote("Linear fit: p < 0.001" ~
                       ", " ~ r^2 ~ " = " ~ .(round(sm1$r.squared, digits = 3)) ~
                       ", AIC = " ~ .(round(AIC(fit1) + (2*3*4/(n-3-1)), digits = 3))))

text(usr[1], usr[3] + 0.1 - 0.25*(usr[3] - usr[4]) / 10, pos = 4, adj = c(1, 1),
     labels = bquote("Quadratic fit:  p < 0.001" ~
                       ", " ~ r^2 ~ " = " ~ .(round(sm2$r.squared, digits = 3)) ~
                       ", AIC = " ~ .(round(AIC(fit2) + (2*4*5/(n-4-1)), digits = 3))))

abline(fit1)


dev.off()

png("figures/qfinal-Fig04d.png", width = 1200, height = 800)
op <- par(mar = c(4, 4, 0.5, 0.5), cex = 2)
#plot.ses(beta.ses$sd1[c(T, F)], mh[c(T, F)], col = "skyblue", lab = expression("SES D'"[pw]), xlab = "Mean altitude, m")
#legend("topright", legend = "Shrub layer", bty = "n")

plot(mh[c(T, F)], beta.ses$sd1[c(T, F)], pch = 19, col = "skyblue", xlab = "Mean altitude, m", ylab = expression("SES D'"[pw]))
fit1 <- lm(beta.ses$sd1[c(T, F)] ~ mh[c(T, F)])
fit2 <- lm(beta.ses$sd1[c(T, F)] ~ mh[c(T, F)] + I(mh[c(T, F)]^2))
sm1 <- summary(fit1)
sm2 <- summary(fit2)
p1 <- pf(sm1$fstatistic[1], df1 = sm1$fstatistic[2], df2 = sm1$fstatistic[3], lower.tail = F)
p2 <- pf(sm2$fstatistic[1], df1 = sm2$fstatistic[2], df2 = sm2$fstatistic[3], lower.tail = F)

usr <- par("usr")

n <- length(beta.ses$sd1[c(T, F)])
text(usr[2], usr[3] + 0.1 - 1.25*(usr[3] - usr[4]) / 10, pos = 2, adj = c(1, 1),
     labels = "Shrub layer")

text(usr[2], usr[3] + 0.1 - 0.75*(usr[3] - usr[4]) / 10, pos = 2, adj = c(1, 1),
     labels = bquote("Linear fit: p = " ~ .(round(p1, digits = 3)) ~
                       ", " ~ r^2 ~ " = " ~ .(round(sm1$r.squared, digits = 3)) ~
                       ", AIC = " ~ .(round(AIC(fit1) + (2*3*4/(n-3-1)), digits = 3))))

text(usr[2], usr[3] + 0.1 - 0.25*(usr[3] - usr[4]) / 10, pos = 2, adj = c(1, 1),
     labels = bquote("Quadratic fit:  p = " ~ .(round(p2, digits = 3)) ~
                       ", " ~ r^2 ~ " = " ~ .(round(sm2$r.squared, digits = 3)) ~
                       ", AIC = " ~ .(round(AIC(fit2) + (2*4*5/(n-4-1)), digits = 3))))


dev.off()

png("figures/qfinal-Fig04e.png", width = 1200, height = 800)
op <- par(mar = c(4, 4, 0.5, 0.5), cex = 2)
#plot.ses(beta.ses$hd0[c(T, F)], mh[c(T, F)], col = "forestgreen", lab = expression("SES " * D[pw]), xlab = "Mean altitude, m")
#legend("topright", legend = "Herb layer", bty = "n")

plot(mh[c(T, F)], beta.ses$hd0[c(T, F)], pch = 19, col = "forestgreen", xlab = "Mean altitude, m", ylab = expression("SES " * D[pw]), ylim = c(-0.22, 2.5))
fit1 <- lm(beta.ses$hd0[c(T, F)] ~ mh[c(T, F)])
fit2 <- lm(beta.ses$hd0[c(T, F)] ~ mh[c(T, F)] + I(mh[c(T, F)]^2))
sm1 <- summary(fit1)
sm2 <- summary(fit2)
p1 <- pf(sm1$fstatistic[1], df1 = sm1$fstatistic[2], df2 = sm1$fstatistic[3], lower.tail = F)
p2 <- pf(sm2$fstatistic[1], df1 = sm2$fstatistic[2], df2 = sm2$fstatistic[3], lower.tail = F)

usr <- par("usr")

n <- length(beta.ses$hd0[c(T, F)])
text(usr[1], usr[4] - 0.25*(usr[4] - usr[3]) / 10, pos = 4, adj = c(1, 1),
     labels = "Herb layer")

text(usr[1], usr[4] - 0.75*(usr[4] - usr[3]) / 10, pos = 4, adj = c(1, 1),
     labels = bquote("Linear fit: p = " ~ .(round(p1, digits = 3)) ~
                       ", " ~ r^2 ~ " = " ~ .(round(sm1$r.squared, digits = 3)) ~
                       ", AIC = " ~ .(round(AIC(fit1) + (2*3*4/(n-3-1)), digits = 3))))

text(usr[1], usr[4] - 1.25*(usr[4] - usr[3]) / 10, pos = 4, adj = c(1, 1),
     labels = bquote("Quadratic fit:  p = " ~ .(round(p2, digits = 3)) ~
                       ", " ~ r^2 ~ " = " ~ .(round(sm2$r.squared, digits = 3)) ~
                       ", AIC = " ~ .(round(AIC(fit2) + (2*4*5/(n-4-1)), digits = 3))))

dev.off()

png("figures/qfinal-Fig04f.png", width = 1200, height = 800)
op <- par(mar = c(4, 4, 0.5, 0.5), cex = 2)
#plot.ses(beta.ses$hd1[c(T, F)], mh[c(T, F)], col = "forestgreen", lab = expression("SES D'"[pw]), xlab = "Mean altitude, m")
#legend("topright", legend = "Herb layer", bty = "n")

plot(mh[c(T, F)], beta.ses$hd1[c(T, F)], pch = 19, col = "forestgreen", xlab = "Mean altitude, m", ylab = expression("SES D'"[pw]), ylim = c(-1.22, 2.7))
fit1 <- lm(beta.ses$hd1[c(T, F)] ~ mh[c(T, F)])
fit2 <- lm(beta.ses$hd1[c(T, F)] ~ mh[c(T, F)] + I(mh[c(T, F)]^2))
sm1 <- summary(fit1)
sm2 <- summary(fit2)
p1 <- pf(sm1$fstatistic[1], df1 = sm1$fstatistic[2], df2 = sm1$fstatistic[3], lower.tail = F)
p2 <- pf(sm2$fstatistic[1], df1 = sm2$fstatistic[2], df2 = sm2$fstatistic[3], lower.tail = F)

usr <- par("usr")

n <- length(beta.ses$hd1[c(T, F)])
text(usr[2], usr[4] - 0.25*(usr[4] - usr[3]) / 10, pos = 2, adj = c(1, 1),
     labels = "Herb layer")

text(usr[2], usr[4] - 0.75*(usr[4] - usr[3]) / 10, pos = 2, adj = c(1, 1),
     labels = bquote("Linear fit: p = " ~ .(round(p1, digits = 3)) ~
                       ", " ~ r^2 ~ " = " ~ .(round(sm1$r.squared, digits = 3)) ~
                       ", AIC = " ~ .(round(AIC(fit1) + (2*3*4/(n-3-1)), digits = 3))))

text(usr[2], usr[4] - 1.25*(usr[4] - usr[3]) / 10, pos = 2, adj = c(1, 1),
     labels = bquote("Quadratic fit:  p = " ~ .(round(p2, digits = 3)) ~
                       ", " ~ r^2 ~ " = " ~ .(round(sm2$r.squared, digits = 3)) ~
                       ", AIC = " ~ .(round(AIC(fit2) + (2*4*5/(n-4-1)), digits = 3))))
x <- seq(min(hgt), max(hgt), len = 1000)
abc <- coef(fit2)
y <- abc[1] + abc[2] * x + abc[3] * x^2
lines(x, y)


dev.off()
