hgt <- read.table("raw.data/hgt.txt")[[1]]

load("clean.data/tree-trait-ses.rda")
load("clean.data/shrub-trait-ses.rda")
load("clean.data/herb-trait-ses.rda")

x <- seq(1050, 1800, len = 1000)

png("figures/report/tree-func-hgt.png", width = 600, height = 400)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot(hgt, t.ses$mpd.z, type = "n", xlab = "Высота, м", ylab = "SES MPD")
abline(h = 0, lwd = 2, col = "grey")
points(hgt, t.ses$mpd.z, pch = 19, col = "tomato")

fit2 <- lm(t.ses$mpd.z ~ hgt + I(hgt^2))
abc <- coef(fit2)
y <- abc[1] + abc[2] * x + abc[3] * x^2
lines(x, y, lwd = 2)
dev.off()

png("figures/report/shrub-func-hgt.png", width = 600, height = 400)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot(hgt, s.ses$mpd.z, type = "n", xlab = "Высота, м", ylab = "SES MPD", ylim = c(-2, 4))
abline(h = 0, lwd = 2, col = "grey")
points(hgt, s.ses$mpd.z, pch = 19, col = "skyblue")

fit2 <- lm(s.ses$mpd.z ~ hgt + I(hgt^2))
abc <- coef(fit2)
y <- abc[1] + abc[2] * x + abc[3] * x^2
lines(x, y, lwd = 2)
dev.off()

png("figures/report/herb-func-hgt.png", width = 600, height = 400)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot(hgt, h.ses$mpd.z, type = "n", xlab = "Высота, м", ylab = "SES MPD")
abline(h = 0, lwd = 2, col = "grey")
points(hgt, h.ses$mpd.z, pch = 19, col = "forestgreen")
dev.off()

png("figures/report/tree-func-hgt-a.png", width = 600, height = 400)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot(hgt, t.ses$mpd.a.z, type = "n", xlab = "Высота, м", ylab = "SES MPD")
abline(h = 0, lwd = 2, col = "grey")
points(hgt, t.ses$mpd.a.z, pch = 19, col = "tomato")
dev.off()

png("figures/report/shrub-func-hgt-a.png", width = 600, height = 400)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot(hgt, s.ses$mpd.a.z, type = "n", xlab = "Высота, м", ylab = "SES MPD")
abline(h = 0, lwd = 2, col = "grey")
points(hgt, s.ses$mpd.a.z, pch = 19, col = "skyblue")

fit1 <- lm(s.ses$mpd.a.z ~ hgt)
abline(fit1, lwd = 2)
dev.off()

png("figures/report/herb-func-hgt-a.png", width = 600, height = 400)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot(hgt, h.ses$mpd.a.z, type = "n", xlab = "Высота, м", ylab = "SES MPD")
abline(h = 0, lwd = 2, col = "grey")
points(hgt, h.ses$mpd.a.z, pch = 19, col = "forestgreen")
dev.off()






load("clean.data/functional-scaling-ses.rda")
sc <- 1:20

png("figures/report/tree-func-scaling.png", width = 600, height = 400)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
meanSES <- sapply(ta_sc_ses, function(x) mean(x[, "z.tl"], na.rm = T))
plot(sc, meanSES, pch = 21, bg = "tomato", xlab = "Уровень агрегации", ylab = "Средний SES MPD")
model <- lm(meanSES ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model, lwd = 2)
dev.off()

png("figures/report/tree-func-scaling-a.png", width = 600, height = 400)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
meanSES <- sapply(ta_sc_ses, function(x) mean(x[, "z.a.tl"], na.rm = T))
plot(sc, meanSES, pch = 21, bg = "tomato", xlab = "Уровень агрегации", ylab = "Средний SES MPD")
model <- lm(meanSES ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model, lwd = 2)
dev.off()



png("figures/report/shrub-func-scaling.png", width = 600, height = 400)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
meanSES <- sapply(sa_sc_ses, function(x) mean(x[, "z.tl"], na.rm = T))
plot(sc, meanSES, pch = 21, bg = "skyblue", xlab = "Уровень агрегации", ylab = "Средний SES MPD")
model <- lm(meanSES ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model, lwd = 2)
dev.off()

png("figures/report/shrub-func-scaling-a.png", width = 600, height = 400)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
meanSES <- sapply(sa_sc_ses, function(x) mean(x[, "z.a.tl"], na.rm = T))
plot(sc, meanSES, pch = 21, bg = "skyblue", xlab = "Уровень агрегации", ylab = "Средний SES MPD")
model <- lm(meanSES ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model, lwd = 2)
dev.off()



png("figures/report/herb-func-scaling.png", width = 600, height = 400)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
meanSES <- sapply(ha_sc_ses, function(x) mean(x[, "z.tl"], na.rm = T))
plot(sc, meanSES, pch = 21, bg = "forestgreen", xlab = "Уровень агрегации", ylab = "Средний SES MPD")
model <- lm(meanSES ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model, lwd = 2)
dev.off()

png("figures/report/herb-func-scaling-a.png", width = 600, height = 400)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
meanSES <- sapply(ha_sc_ses, function(x) mean(x[, "z.a.tl"], na.rm = T))
plot(sc, meanSES, pch = 21, bg = "forestgreen", xlab = "Уровень агрегации", ylab = "Средний SES MPD")
model <- lm(meanSES ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model, lwd = 2)
dev.off()