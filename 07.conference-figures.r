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

ta <- ta[, order(colSums(ta))]
sa <- sa[, order(colSums(sa))]
ha <- ha[, order(colSums(ha))]


t.col <- ta
t.col[] <- "white"
t.col[ta > 0 & ta < 15] <- "tomato"
t.col[ta >= 15 & ta < 30] <- "darkred"
t.col[ta >= 30] <- "black"

ab <- colSums(ta)
ab <- ab / max(ab) * 40


png("figures/conf-trees.png", width = 1200, height = 600, bg = "transparent")

op <- par(mar = c(0.25, 10, 0.25, 0.25), cex = 2)

plot(1, 1, type = "n", ylim = c(-1, 16.5), xlim = c(0, 151), axes = F, ann = F)
for (ii in 1:16) points(1:96, rep(ii, 96), pch = 22, cex = 1, bg = t.col[,ii])

rect(xleft = 111, ybottom = (1:16) - 0.2, xright = 111 + ab, ytop = (1:16) + 0.2, col = "tomato")

axis(side = 2, at = c(1:16), labels = sub("_" , " ", colnames(ta)), las = 2, tick = F)
legend("bottom", c("0", "1-15", "16-30", "> 30"), horiz = TRUE, bty = "n", pch = 22, 
       pt.bg = c("white", "tomato", "darkred", "black"), pt.cex = 2)

dev.off()

max(sa)

s.col <- sa
s.col[] <- "white"
s.col[sa > 0 & sa < 10] <- "skyblue"
s.col[sa >= 10 & sa < 100] <- "darkblue"
s.col[sa >= 100] <- "black"

ab <- colSums(sa)
ab <- ab / max(ab) * 40

png("figures/conf-shrubs.png", width = 1200, height = 600, bg = "transparent")

op <- par(mar = c(0.25, 12, 0.25, 0.25), cex = 2)

plot(1, 1, type = "n", ylim = c(-1, 16.5), xlim = c(0, 151), axes = F, ann = F)
for (ii in 1:16) points(1:96, rep(ii, 96), pch = 22, cex = 1, bg = s.col[,ii+19])

rect(xleft = 111, ybottom = (1:16) - 0.2, xright = 111 + ab[20:35], ytop = (1:16) + 0.2, col = "skyblue")

axis(side = 2, at = c(1:16), labels = sub("_" , " ", colnames(sa)[20:35]), las = 2, tick = F)
legend("bottom", c("0", "1-10", "10-100", "> 100"), horiz = TRUE, bty = "n", pch = 22, 
       pt.bg = c("white", "skyblue", "darkblue", "black"), pt.cex = 2)

dev.off()





max(ha)

h.col <- ha
h.col[] <- "white"
h.col[ha > 0 & ha < 10] <- "lightgreen"
h.col[ha >= 10 & ha < 50] <- "darkgreen"
h.col[ha >= 50] <- "black"

ab <- colSums(ha)
ab <- ab / max(ab) * 40

png("figures/conf-herbs.png", width = 1200, height = 600, bg = "transparent")

op <- par(mar = c(0.25, 12, 0.25, 0.25), cex = 2)

plot(1, 1, type = "n", ylim = c(-1, 16.5), xlim = c(0, 151), axes = F, ann = F)
for (ii in 1:16) points(1:96, rep(ii, 96), pch = 22, cex = 1, bg = h.col[,ii+133])

rect(xleft = 111, ybottom = (1:16) - 0.2, xright = 111 + ab[134:149], ytop = (1:16) + 0.2, col = "lightgreen")

axis(side = 2, at = c(1:16), labels = sub("_" , " ", colnames(ha)[134:149]), las = 2, tick = F)
legend("bottom", c("0", "1-10", "10-100", "> 100"), horiz = TRUE, bty = "n", pch = 22, 
       pt.bg = c("white", "lightgreen", "darkgreen", "black"), pt.cex = 2)

dev.off()

#---#

x <- seq(1050, 1800, len = 1000)

png("figures/conf-2a-tree-cov.png", width = 800, height = 550)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot(hgt, rowSums(ta), pch = 19, col = "tomato", 
     xlab = "Высота, м", ylab = "N, число стволов")
fit2 <- lm(rowSums(ta) ~ hgt + I(hgt^2))
abc <- coef(fit2)
y <- abc[1] + abc[2] * x + abc[3] * x^2
lines(x, y, lwd = 2)
legend("topright", legend = "Деревья", bty = "n")
dev.off()

png("figures/conf-2b-shrub-cov.png", width = 800, height = 550)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot(hgt, rowSums(sa), pch = 19, col = "skyblue", 
     xlab = "Высота, м", ylab = "N, число стволов")
fit2 <- lm(rowSums(sa) ~ hgt + I(hgt^2))
abc <- coef(fit2)
y <- abc[1] + abc[2] * x + abc[3] * x^2
lines(x, y, lwd = 2)
legend("topright", legend = "Кустарники", bty = "n")
dev.off()

png("figures/conf-2с-herb-cov.png", width = 800, height = 550)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot(hgt, rowSums(ha), pch = 19, col = "forestgreen", 
     xlab = "Высота, м", ylab = "N, проект. покр.")
fit2 <- lm(rowSums(ha) ~ hgt + I(hgt^2))
abc <- coef(fit2)
y <- abc[1] + abc[2] * x + abc[3] * x^2
lines(x, y, lwd = 2)
legend("topright", legend = "Трава", bty = "n")
dev.off()

png("figures/conf-3a-tree-S.png", width = 800, height = 550)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot(hgt, specnumber(ta), pch = 19, col = "tomato", 
     xlab = "Высота, м", ylab = "S, число видов")
fit2 <- lm(specnumber(ta) ~ hgt + I(hgt^2))
abc <- coef(fit2)
y <- abc[1] + abc[2] * x + abc[3] * x^2
lines(x, y, lwd = 2)
legend("topright", legend = "Деревья", bty = "n")
dev.off()

png("figures/conf-3b-shrub-S.png", width = 800, height = 550)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot(hgt, specnumber(sa), pch = 19, col = "skyblue", 
     xlab = "Высота, м", ylab = "S, число видов")
fit1 <- lm(specnumber(sa) ~ hgt)
abline(fit1, lwd = 2)
legend("topright", legend = "Кустарники", bty = "n")
dev.off()

png("figures/conf-3c-herb-S.png", width = 800, height = 550)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot(hgt, specnumber(ha), pch = 19, col = "forestgreen", 
     xlab = "Высота, м", ylab = "S, число видов")
fit2 <- lm(specnumber(ha) ~ hgt + I(hgt^2))
abc <- coef(fit2)
y <- abc[1] + abc[2] * x + abc[3] * x^2
lines(x, y, lwd = 2)
legend("topleft", legend = "Трава", bty = "n")
dev.off()


png("figures/conf-4a-tree-phylo.png", width = 800, height = 600)
par(mar = c(2, 0, 0, 0))
plot.phylo(wood.tree)
axisPhylo()
dev.off()

png("figures/conf-4b-shrub-phylo.png", width = 800, height = 600)
par(mar = c(2, 0, 0, 0))
plot.phylo(shrub.tree)
axisPhylo()
dev.off()

png("figures/conf-4c-herb-phylo.png", width = 800, height = 600)
par(mar = c(2, 0, 0, 0))
plot.phylo(herb.tree, cex = 0.5)
axisPhylo()
dev.off()

library(ecomf)
tpd <- PqD(ta, wood.tree, hill = F)
spd <- PqD(sa, shrub.tree, hill = F)
hpd <- PqD(ha, herb.tree, hill = F)

png("figures/conf-5a-tree-PD0.png", width = 700, height = 400)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot(hgt, tpd[, 1], pch = 19, col = "tomato", 
     xlab = "Высота, м", ylab = expression(PD[0]))
fit2 <- lm(tpd[, 1] ~ hgt + I(hgt^2))
abc <- coef(fit2)
y <- abc[1] + abc[2] * x + abc[3] * x^2
lines(x, y, lwd = 2)
dev.off()

png("figures/conf-5b-tree-PD2.png", width = 700, height = 400)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot(hgt, tpd[, 3], pch = 19, col = "tomato", 
     xlab = "Высота, м", ylab = expression(PD[2]))
fit2 <- lm(tpd[, 3] ~ hgt + I(hgt^2))
abc <- coef(fit2)
y <- abc[1] + abc[2] * x + abc[3] * x^2
lines(x, y, lwd = 2)
dev.off()

png("figures/conf-5c-shrub-PD0.png", width = 700, height = 400)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot(hgt, spd[, 1], pch = 19, col = "skyblue", 
     xlab = "Высота, м", ylab = expression(PD[0]))
fit1 <- lm(spd[, 1] ~ hgt)
abline(fit1, lwd = 2)
dev.off()

png("figures/conf-5d-shrub-PD2.png", width = 700, height = 400)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot(hgt, spd[, 3], pch = 19, col = "skyblue", 
     xlab = "Высота, м", ylab = expression(PD[2]))
fit1 <- lm(spd[, 3] ~ hgt)
abline(fit1, lwd = 2)
dev.off()

png("figures/conf-5e-herb-PD0.png", width = 700, height = 400)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot(hgt, hpd[, 1], pch = 19, col = "forestgreen", 
     xlab = "Высота, м", ylab = expression(PD[0]))
fit2 <- lm(hpd[, 1] ~ hgt + I(hgt^2))
abc <- coef(fit2)
y <- abc[1] + abc[2] * x + abc[3] * x^2
lines(x, y, lwd = 2)
dev.off()

png("figures/conf-5f-herb-PD2.png", width = 700, height = 400)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot(hgt, hpd[, 3], pch = 19, col = "forestgreen", 
     xlab = "Высота, м", ylab = expression(PD[2]))
fit2 <- lm(hpd[, 3] ~ hgt + I(hgt^2))
abc <- coef(fit2)
y <- abc[1] + abc[2] * x + abc[3] * x^2
lines(x, y, lwd = 2)
dev.off()


#---#

null <- replicate(999, mpd(randomizeMatrix(ha, null.model = "independentswap"), herb.mat))
nri <- -(mpd(ha, herb.mat) - rowMeans(null)) / apply(null, 1, sd)

plot.ses(nri, hgt, "tomato", "nri")

which.min(nri)

hist(null[which.min(nri), ])

hist(null[which.max(nri), ])

png("figures/conf-6a-herb-hist1.png", width = 600, height = 400)
par(mar = c(4,4,.5,.5), cex = 2)
num <- which.min(nri)
obs <- mpd(ha[num, ], herb.mat)
h <- hist(null[num, ], breaks = 50, col = "limegreen", main = "", xlab = "MPD", ylab = "Частота")
arrows(x0 = obs, y0 = max(h$counts)/1.5, y1 = 30, angle = 25, lwd = 2)
text(obs-10, max(h$counts)/1.5, labels = bquote(MPD[obs] == .(round(obs, digits = 1))), pos = 3 )
dev.off()

mean(null[num, ])
sd(null[num, ])
hgt[num]

png("figures/conf-6b-herb-hist2.png", width = 600, height = 400)
par(mar = c(4,4,.5,.5), cex = 2)
num <- which.max(nri)
obs <- mpd(ha[num, ], herb.mat)
h <- hist(null[num, ], breaks = 50, col = "limegreen", main = "", xlab = "MPD", ylab = "Частота")
arrows(x0 = obs, y0 = max(h$counts)/1.15, y1 = 35, angle = 25, lwd = 2)
text(obs, max(h$counts)/1.15, labels = bquote(MPD[obs] == .(round(obs, digits = 1))), pos = 3 )
dev.off()

mean(null[num, ])
sd(null[num, ])
hgt[num]


#---#

load("clean.data/phylo-ses-pqd.rda")

png("figures/conf-7a-tree-NRI.png", width = 800, height = 550)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot(hgt, -ta.ses$mpd, type = "n", xlab = "Высота, м", ylab = expression(NRI))
abline(h = 0, lwd = 2, col = "grey")
points(hgt, -ta.ses$mpd, pch = 19, col = "tomato")
fit1 <- lm(-ta.ses$mpd ~ hgt)
abline(fit1, lwd = 2)
legend("topright", legend = "Деревья", bty = "n")
dev.off()

png("figures/conf-7b-shrub-NRI.png", width = 800, height = 550)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot(hgt, -sa.ses$mpd, type = "n", xlab = "Высота, м", ylab = expression(NRI), ylim = c(-2, 4))
abline(h = 0, lwd = 2, col = "grey")
points(hgt, -sa.ses$mpd, pch = 19, col = "skyblue")
fit1 <- lm(-sa.ses$mpd ~ hgt)
abline(fit1, lwd = 2)
legend("topleft", legend = "Кустарники", bty = "n")
dev.off()

png("figures/conf-7c-herb-NRI.png", width = 800, height = 550)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot(hgt, -ha.ses$mpd, type = "n", xlab = "Высота, м", ylab = expression(NRI))
abline(h = 0, lwd = 2, col = "grey")
points(hgt, -ha.ses$mpd, pch = 19, col = "forestgreen")
fit1 <- lm(-ha.ses$mpd ~ hgt)
abline(fit1, lwd = 2)
legend("bottomright", legend = "Трава", bty = "n")
dev.off()

png("figures/conf-8a-tree-NRIa.png", width = 800, height = 550)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot(hgt, -ta.ses$mpd.a, type = "n", xlab = "Высота, м", ylab = expression(NRI[a]))
abline(h = 0, lwd = 2, col = "grey")
points(hgt, -ta.ses$mpd.a, pch = 19, col = "tomato")
fit1 <- lm(-ta.ses$mpd.a ~ hgt)
abline(fit1, lwd = 2)
dev.off()

png("figures/conf-8b-shrub-NRIa.png", width = 800, height = 550)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot(hgt, -sa.ses$mpd.a, type = "n", xlab = "Высота, м", ylab = expression(NRI[a]))
abline(h = 0, lwd = 2, col = "grey")
points(hgt, -sa.ses$mpd.a, pch = 19, col = "skyblue")
fit2 <- lm(-sa.ses$mpd.a ~ hgt + I(hgt^2))
abc <- coef(fit2)
y <- abc[1] + abc[2] * x + abc[3] * x^2
lines(x, y, lwd = 2)
dev.off()

png("figures/conf-8c-herb-NRIa.png", width = 800, height = 550)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot(hgt, -ha.ses$mpd.a, type = "n", xlab = "Высота, м", ylab = expression(NRI[a]))
abline(h = 0, lwd = 2, col = "grey")
points(hgt, -ha.ses$mpd.a, pch = 19, col = "forestgreen")
fit2 <- lm(-ha.ses$mpd.a ~ hgt + I(hgt^2))
abc <- coef(fit2)
y <- abc[1] + abc[2] * x + abc[3] * x^2
lines(x, y, lwd = 2)
dev.off()



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

#---#

load("clean.data/phylo-ses-pqd.rda")

ses <- ha.ses.sc[[1]]$pqd.1
lab <- "PqD"
col <- "forestgreen"

plot.ses <- function(ses, hgt = hgt, col, lab) {
  plot(hgt, ses, pch = 19, col = col, ylab = lab)
  fit1 <- lm(ses ~ hgt)
  fit2 <- lm(ses ~ hgt + I(hgt^2))
  sm1 <- summary(fit1)
  sm2 <- summary(fit2)
  p1 <- pf(sm1$fstatistic[1], df1 = sm1$fstatistic[2], df2 = sm1$fstatistic[3], lower.tail = F)
  p2 <- pf(sm2$fstatistic[1], df1 = sm2$fstatistic[2], df2 = sm2$fstatistic[3], lower.tail = F)
  usr <- par("usr")
  if (p1 < 0.05)
  {
    if (p2 >= 0.05)
    {
      abline(fit1)
      text(usr[1], usr[3] + (usr[4] - usr[3]) / 10, round(sm1$r.squared, digits = 3), pos = 4, adj = c(1, 1))
    }
    else
    {
      n <- length(ses)
      delta <- AIC(fit2) + (2*4*5/(n-4-1)) - AIC(fit1) - (2*3*4/(n-3-1))
      if (delta < 0)
      {
        x <- seq(min(hgt), max(hgt), len = 1000)
        abc <- coef(fit2)
        y <- abc[1] + abc[2] * x + abc[3] * x^2
        lines(x, y)
        text(usr[1], usr[3] + (usr[4] - usr[3]) / 10, round(sm2$r.squared, digits = 3), pos = 4, adj = c(1, 1))
      }
      else
      {
        abline(fit1)
        text(usr[1], usr[3] + (usr[4] - usr[3]) / 10, round(sm1$r.squared, digits = 3), pos = 4, adj = c(1, 1))
      }
    }
  }
  else if (p2 < 0.05)
  {
    x <- seq(min(hgt), max(hgt), len = 1000)
    abc <- coef(fit2)
    y <- abc[1] + abc[2] * x + abc[3] * x^2
    lines(x, y)
    text(usr[1], usr[3] + (usr[4] - usr[3]) / 10, round(sm2$r.squared, digits = 3), pos = 4, adj = c(1, 1))
  }
    
}

op <- par(mfcol = c(7, 3), mar = c(4,4,0.5,0.5))
plot.ses(ta.ses$mpd, hgt, col = "tomato", lab = "MPD")
plot.ses(ta.ses$mpd.a, hgt, col = "tomato", lab = "MPD.a")
plot.ses(ta.ses$pqd.n1, hgt, col = "tomato", lab = expression(PD[-1]))
plot.ses(ta.ses$pqd.0, hgt, col = "tomato", lab = expression(PD[0]))
plot.ses(ta.ses$pqd.1, hgt, col = "tomato", lab = expression(PD[1]))
plot.ses(ta.ses$pqd.2, hgt, col = "tomato", lab = expression(PD[2]))
plot.ses(ta.ses$pqd.5, hgt, col = "tomato", lab = expression(PD[5]))

plot.ses(sa.ses$mpd, hgt, col = "skyblue", lab = "MPD")
plot.ses(sa.ses$mpd.a, hgt, col = "skyblue", lab = "MPD.a")
plot.ses(sa.ses$pqd.n1, hgt, col = "skyblue", lab = expression(PD[-1]))
plot.ses(sa.ses$pqd.0, hgt, col = "skyblue", lab = expression(PD[0]))
plot.ses(sa.ses$pqd.1, hgt, col = "skyblue", lab = expression(PD[1]))
plot.ses(sa.ses$pqd.2, hgt, col = "skyblue", lab = expression(PD[2]))
plot.ses(sa.ses$pqd.5, hgt, col = "skyblue", lab = expression(PD[5]))

plot.ses(ha.ses$mpd, hgt, col = "forestgreen", lab = "MPD")
plot.ses(ha.ses$mpd.a, hgt, col = "forestgreen", lab = "MPD.a")
plot.ses(ha.ses$pqd.n1, hgt, col = "forestgreen", lab = expression(PD[-1]))
plot.ses(ha.ses$pqd.0, hgt, col = "forestgreen", lab = expression(PD[0]))
plot.ses(ha.ses$pqd.1, hgt, col = "forestgreen", lab = expression(PD[1]))
plot.ses(ha.ses$pqd.2, hgt, col = "forestgreen", lab = expression(PD[2]))
plot.ses(ha.ses$pqd.5, hgt, col = "forestgreen", lab = expression(PD[5]))
par(op)

op <- par(mfrow = c(3, 2), mar = c(4,4,0.5,0.5))
plot.ses(ta.ses$mpd, col = "tomato", lab = "MPD")
plot.ses(ta.ses.sc[[1]]$mpd, col = "tomato", lab = "MPD")
plot.ses(sa.ses$mpd, col = "skyblue", lab = "MPD")
plot.ses(sa.ses.sc[[1]]$mpd, col = "skyblue", lab = "MPD")
plot.ses(ha.ses$mpd, col = "forestgreen", lab = "MPD")
plot.ses(ha.ses.sc[[1]]$mpd, col = "forestgreen", lab = "MPD")
par(op)

op <- par(mfrow = c(3, 2), mar = c(4,4,0.5,0.5))
plot(ta.ses$mpd, ta.ses$pqd.0, pch = 19, col = "tomato")
abline(0, 1)
plot(ta.ses$mpd.a, ta.ses$pqd.1, pch = 19, col = "tomato")
abline(0, 1)
plot(sa.ses$mpd, sa.ses$pqd.0, pch = 19, col = "skyblue")
abline(0, 1)
plot(sa.ses$mpd.a, sa.ses$pqd.1, pch = 19, col = "skyblue")
abline(0, 1)
plot(ha.ses$mpd, ha.ses$pqd.0, pch = 19, col = "forestgreen")
abline(0, 1)
plot(ha.ses$mpd.a, ha.ses$pqd.1, pch = 19, col = "forestgreen")
abline(0, 1)
par(op)


sc <- 1:20

op <- par(mfcol = c(7, 3), mar = c(4,4,0.5,0.5))
plot.ses(sapply(ta.ses.sc, function(x) mean(x$mpd, na.rm = T)), sc, col = "tomato", lab = "MPD")
plot.ses(sapply(ta.ses.sc, function(x) mean(x$mpd.a, na.rm = T)), sc, col = "tomato", lab = "MPD")
plot.ses(sapply(ta.ses.sc, function(x) mean(x$pqd.n1, na.rm = T)), sc, col = "tomato", lab = "PD-1")
plot.ses(sapply(ta.ses.sc, function(x) mean(x$pqd.0, na.rm = T)), sc, col = "tomato", lab = "PD0")
plot.ses(sapply(ta.ses.sc, function(x) mean(x$pqd.1, na.rm = T)), sc, col = "tomato", lab = "PD1")
plot.ses(sapply(ta.ses.sc, function(x) mean(x$pqd.2, na.rm = T)), sc, col = "tomato", lab = "PD2")
plot.ses(sapply(ta.ses.sc, function(x) mean(x$pqd.5, na.rm = T)), sc, col = "tomato", lab = "PD5")

plot.ses(sapply(sa.ses.sc, function(x) mean(x$mpd, na.rm = T)), sc, col = "skyblue", lab = "MPD")
plot.ses(sapply(sa.ses.sc, function(x) mean(x$mpd.a, na.rm = T)), sc, col = "skyblue", lab = "MPD")
plot.ses(sapply(sa.ses.sc, function(x) mean(x$pqd.n1, na.rm = T)), sc, col = "skyblue", lab = "PD-1")
plot.ses(sapply(sa.ses.sc, function(x) mean(x$pqd.0, na.rm = T)), sc, col = "skyblue", lab = "PD0")
plot.ses(sapply(sa.ses.sc, function(x) mean(x$pqd.1, na.rm = T)), sc, col = "skyblue", lab = "PD1")
plot.ses(sapply(sa.ses.sc, function(x) mean(x$pqd.2, na.rm = T)), sc, col = "skyblue", lab = "PD2")
plot.ses(sapply(sa.ses.sc, function(x) mean(x$pqd.5, na.rm = T)), sc, col = "skyblue", lab = "PD5")

plot.ses(sapply(ha.ses.sc, function(x) mean(x$mpd, na.rm = T)), sc, col = "forestgreen", lab = "MPD")
plot.ses(sapply(ha.ses.sc, function(x) mean(x$mpd.a, na.rm = T)), sc, col = "forestgreen", lab = "MPD")
plot.ses(sapply(ha.ses.sc, function(x) mean(x$pqd.n1, na.rm = T)), sc, col = "forestgreen", lab = "PD-1")
plot.ses(sapply(ha.ses.sc, function(x) mean(x$pqd.0, na.rm = T)), sc, col = "forestgreen", lab = "PD0")
plot.ses(sapply(ha.ses.sc, function(x) mean(x$pqd.1, na.rm = T)), sc, col = "forestgreen", lab = "PD1")
plot.ses(sapply(ha.ses.sc, function(x) mean(x$pqd.2, na.rm = T)), sc, col = "forestgreen", lab = "PD2")
plot.ses(sapply(ha.ses.sc, function(x) mean(x$pqd.5, na.rm = T)), sc, col = "forestgreen", lab = "PD5")
par(op)


plot(sapply(ta.ses.sc, function(x) mean(x$mpd.a, na.rm = T)), 
     sapply(ta.ses.sc, function(x) mean(x$pqd.1, na.rm = T)), pch = 19, col = "tomato")
abline(0, 1)

plot(sapply(sa.ses.sc, function(x) mean(x$mpd.a, na.rm = T)), 
     sapply(sa.ses.sc, function(x) mean(x$pqd.1, na.rm = T)), pch = 19, col = "skyblue")
abline(0, 1)

plot(-sapply(ha.ses.sc, function(x) mean(x$mpd.a, na.rm = T)), 
     sapply(ha.ses.sc, function(x) mean(x$pqd.1, na.rm = T)), pch = 19, col = "forestgreen")
abline(0, 1)


#------------------------------------------------------------------------#

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


tpd <- PqD(ta, wood.tree, q = -1:2, hill = F)
spd <- PqD(sa, shrub.tree, q = -1:2, hill = F)
hpd <- PqD(ha, herb.tree, q = -1:2, hill = F)

op <- par(mfcol = c(7, 3), mar = c(4,4,0.5,0.5))
plot.ses(specnumber(ta), hgt, col = "tomato", lab = expression(D[0] == S))
plot.ses(exp(diversity(ta)), hgt, col = "tomato", lab = expression(D[1]))
plot.ses(diversity(ta, "invsimpson"), hgt, col = "tomato", lab = expression(D[2]))
plot.ses(tpd[, 1], hgt, col = "tomato", lab = expression(PD[-1]))
plot.ses(tpd[, 2], hgt, col = "tomato", lab = expression(PD[0]))
plot.ses(tpd[, 3], hgt, col = "tomato", lab = expression(PD[1]))
plot.ses(tpd[, 4], hgt, col = "tomato", lab = expression(PD[2]))

plot.ses(specnumber(sa), hgt, col = "skyblue", lab = expression(D[0] == S))
plot.ses(exp(diversity(sa)), hgt, col = "skyblue", lab = expression(D[1]))
plot.ses(diversity(sa, "invsimpson"), hgt, col = "skyblue", lab = expression(D[2]))
plot.ses(spd[, 1], hgt, col = "skyblue", lab = expression(PD[-1]))
plot.ses(spd[, 2], hgt, col = "skyblue", lab = expression(PD[0]))
plot.ses(spd[, 3], hgt, col = "skyblue", lab = expression(PD[1]))
plot.ses(spd[, 4], hgt, col = "skyblue", lab = expression(PD[2]))

plot.ses(specnumber(ha), hgt, col = "forestgreen", lab = expression(D[0] == S))
plot.ses(exp(diversity(ha)), hgt, col = "forestgreen", lab = expression(D[1]))
plot.ses(diversity(ha, "invsimpson"), hgt, col = "forestgreen", lab = expression(D[2]))
plot.ses(hpd[, 1], hgt, col = "forestgreen", lab = expression(PD[-1]))
plot.ses(hpd[, 2], hgt, col = "forestgreen", lab = expression(PD[0]))
plot.ses(hpd[, 3], hgt, col = "forestgreen", lab = expression(PD[1]))
plot.ses(hpd[, 4], hgt, col = "forestgreen", lab = expression(PD[2]))
par(op)



plot.ses(rowSums(ta), hgt, col = "tomato", lab = expression(D[0] == S))
plot.ses(rowSums(sa), hgt, col = "skyblue", lab = expression(D[0] == S))
plot.ses(rowSums(ha), hgt, col = "forestgreen", lab = expression(D[0] == S))
