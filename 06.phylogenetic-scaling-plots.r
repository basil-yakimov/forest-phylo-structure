hgt <- read.table("raw.data/hgt.txt")[[1]]

load("clean.data/phylo-ses-pqd.rda")

x <- seq(1050, 1800, len = 1000)

png("figures/report/tree-phylo-hgt.png", width = 600, height = 400)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot(hgt, ta.ses$mpd, type = "n", xlab = "Высота, м", ylab = "SES MPD")
abline(h = 0, lwd = 2, col = "grey")
points(hgt, ta.ses$mpd, pch = 19, col = "tomato")
fit1 <- lm(ta.ses$mpd ~ hgt)
abline(fit1, lwd = 2)
dev.off()

png("figures/report/shrub-phylo-hgt.png", width = 600, height = 400)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot(hgt, sa.ses$mpd, type = "n", xlab = "Высота, м", ylab = "SES MPD", ylim = c(-2, 4))
abline(h = 0, lwd = 2, col = "grey")
points(hgt, sa.ses$mpd, pch = 19, col = "skyblue")
fit1 <- lm(sa.ses$mpd ~ hgt)
abline(fit1, lwd = 2)
dev.off()

png("figures/report/herb-phylo-hgt.png", width = 600, height = 400)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot(hgt, ha.ses$mpd, type = "n", xlab = "Высота, м", ylab = "SES MPD")
abline(h = 0, lwd = 2, col = "grey")
points(hgt, ha.ses$mpd, pch = 19, col = "forestgreen")
fit1 <- lm(ha.ses$mpd ~ hgt)
abline(fit1, lwd = 2)
dev.off()

png("figures/report/tree-phylo-hgt-a.png", width = 600, height = 400)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot(hgt, ta.ses$mpd.a, type = "n", xlab = "Высота, м", ylab = "SES MPD")
abline(h = 0, lwd = 2, col = "grey")
points(hgt, ta.ses$mpd.a, pch = 19, col = "tomato")
fit1 <- lm(ta.ses$mpd.a ~ hgt)
abline(fit1, lwd = 2)
dev.off()

png("figures/report/shrub-phylo-hgt-a.png", width = 600, height = 400)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot(hgt, sa.ses$mpd.a, type = "n", xlab = "Высота, м", ylab = "SES MPD")
abline(h = 0, lwd = 2, col = "grey")
points(hgt, sa.ses$mpd.a, pch = 19, col = "skyblue")
fit2 <- lm(sa.ses$mpd.a ~ hgt + I(hgt^2))
abc <- coef(fit2)
y <- abc[1] + abc[2] * x + abc[3] * x^2
lines(x, y, lwd = 2)
dev.off()

png("figures/report/herb-phylo-hgt-a.png", width = 600, height = 400)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot(hgt, ha.ses$mpd.a, type = "n", xlab = "Высота, м", ylab = "SES MPD")
abline(h = 0, lwd = 2, col = "grey")
points(hgt, ha.ses$mpd.a, pch = 19, col = "forestgreen")
fit2 <- lm(ha.ses$mpd.a ~ hgt + I(hgt^2))
abc <- coef(fit2)
y <- abc[1] + abc[2] * x + abc[3] * x^2
lines(x, y, lwd = 2)
dev.off()







load("clean.data/scaling.rda")
load("clean.data/scaling-ses.rda")
load("clean.data/source-pool-scaling.rda")
load("clean.data/coherent-scaling.rda")

library(effects)
library(nlme)

res$id <- factor(res$id)


plot.ef <- function(fit, col = "limegreen")
{
  ef <- effect("r", fit)
  plot(ef$x$r, ef$fit, type = "l", lwd = 2, col = col, ylim = range(ef$lower, ef$upper))
  lines(ef$x$r, ef$lower, col = col)
  lines(ef$x$r, ef$upper, col = col)
}


sc <- 1:20
rmax <- floor((96-1)/2)




load("clean.data/scaling-ses.rda")

png("figures/report/tree-phylo-scaling.png", width = 600, height = 400)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
meanSES <- sapply(ta_sc_ses, function(x) mean(x[, "z.tl"], na.rm = T))
plot(sc, meanSES, pch = 21, bg = "tomato", xlab = "Уровень агрегации", ylab = "Средний SES MPD")
model <- lm(meanSES ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model, lwd = 2)
dev.off()

png("figures/report/tree-phylo-scaling-a.png", width = 600, height = 400)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
meanSES <- sapply(ta_sc_ses, function(x) mean(x[, "z.a.tl"], na.rm = T))
plot(sc, meanSES, pch = 21, bg = "tomato", xlab = "Уровень агрегации", ylab = "Средний SES MPD")
model <- lm(meanSES ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model, lwd = 2)
dev.off()


load("clean.data/source-pool-scaling.rda")

png("figures/report/tree-phylo-sp-scaling.png", width = 600, height = 400)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
fit <- lme(fixed = t.ses ~ r, random = ~ 1 | id, data = res, method = "REML", na.action = "na.omit")
ef <- effect("r", fit)
plot(t.ses ~ r, res, col = "tomato", pch = 16, xlab = "Размер видового пула", ylab = "SES MPD")
lines(ef$x$r, ef$fit, type = "l", lwd = 3, ylim = range(ef$lower, ef$upper), pch = 16)
lines(ef$x$r, ef$lower, lwd = 2, col = "grey40")
lines(ef$x$r, ef$upper, lwd = 2, col = "grey40")
dev.off()

png("figures/report/tree-phylo-sp-scaling-a.png", width = 600, height = 400)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
fit <- lme(fixed = t.ses.a ~ r, random = ~ 1 | id, data = res, method = "REML", na.action = "na.omit")
ef <- effect("r", fit)
plot(t.ses.a ~ r, res, col = "tomato", pch = 16, xlab = "Размер видового пула", ylab = "SES MPD")
lines(ef$x$r, ef$fit, type = "l", lwd = 3, ylim = range(ef$lower, ef$upper), pch = 16)
lines(ef$x$r, ef$lower, lwd = 2, col = "grey40")
lines(ef$x$r, ef$upper, lwd = 2, col = "grey40")
dev.off()


load("clean.data/coherent-scaling.rda")

png("figures/report/tree-phylo-coherent-scaling.png", width = 600, height = 400)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot(t.ses ~ sc, res, col = "tomato", pch = 16, xlab = "Уровень агрегации", ylab = "SES MPD")
if (anova(lm(t.ses ~ sc, res))[1, 5] < 0.05) abline(lm(t.ses ~ sc, res), lwd = 2)
dev.off()

png("figures/report/tree-phylo-coherent-scaling-a.png", width = 600, height = 400)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot(t.ses.a ~ sc, res, col = "tomato", pch = 16, xlab = "Уровень агрегации", ylab = "SES MPD")
if (anova(lm(t.ses.a ~ sc, res))[1, 5] < 0.05) abline(lm(t.ses.a ~ sc, res), lwd = 2)
dev.off()





load("clean.data/scaling-ses.rda")

png("figures/report/shrub-phylo-scaling.png", width = 600, height = 400)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
meanSES <- sapply(sa_sc_ses, function(x) mean(x[, "z.tl"], na.rm = T))
plot(sc, meanSES, pch = 21, bg = "skyblue", xlab = "Уровень агрегации", ylab = "Средний SES MPD")
model <- lm(meanSES ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model, lwd = 2)
dev.off()

png("figures/report/shrub-phylo-scaling-a.png", width = 600, height = 400)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
meanSES <- sapply(sa_sc_ses, function(x) mean(x[, "z.a.tl"], na.rm = T))
plot(sc, meanSES, pch = 21, bg = "skyblue", xlab = "Уровень агрегации", ylab = "Средний SES MPD")
model <- lm(meanSES ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model, lwd = 2)
dev.off()

load("clean.data/source-pool-scaling.rda")

png("figures/report/shrub-phylo-sp-scaling.png", width = 600, height = 400)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
fit <- lme(fixed = s.ses ~ r, random = ~ 1 | id, data = res, method = "REML", na.action = "na.omit")
ef <- effect("r", fit)
plot(s.ses ~ r, res, col = "skyblue", pch = 16, xlab = "Размер видового пула", ylab = "SES MPD")
lines(ef$x$r, ef$fit, type = "l", lwd = 3, ylim = range(ef$lower, ef$upper), pch = 16)
lines(ef$x$r, ef$lower, lwd = 2, col = "grey40")
lines(ef$x$r, ef$upper, lwd = 2, col = "grey40")
dev.off()

png("figures/report/shrub-phylo-sp-scaling-a.png", width = 600, height = 400)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
fit <- lme(fixed = s.ses.a ~ r, random = ~ 1 | id, data = res, method = "REML", na.action = "na.omit")
ef <- effect("r", fit)
plot(s.ses.a ~ r, res, col = "skyblue", pch = 16, xlab = "Размер видового пула", ylab = "SES MPD")
lines(ef$x$r, ef$fit, type = "l", lwd = 3, ylim = range(ef$lower, ef$upper), pch = 16)
lines(ef$x$r, ef$lower, lwd = 2, col = "grey40")
lines(ef$x$r, ef$upper, lwd = 2, col = "grey40")
dev.off()

load("clean.data/coherent-scaling.rda")

png("figures/report/shrub-phylo-coherent-scaling.png", width = 600, height = 400)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot(s.ses ~ sc, res, col = "skyblue", pch = 16, xlab = "Уровень агрегации", ylab = "SES MPD")
if (anova(lm(s.ses ~ sc, res))[1, 5] < 0.05) abline(lm(s.ses ~ sc, res), lwd = 2)
dev.off()

png("figures/report/shrub-phylo-coherent-scaling-a.png", width = 600, height = 400)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot(s.ses.a ~ sc, res, col = "skyblue", pch = 16, xlab = "Уровень агрегации", ylab = "SES MPD")
if (anova(lm(s.ses.a ~ sc, res))[1, 5] < 0.05) abline(lm(s.ses.a ~ sc, res), lwd = 2)
dev.off()





load("clean.data/scaling-ses.rda")

png("figures/report/herb-phylo-scaling.png", width = 600, height = 400)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
meanSES <- sapply(ha_sc_ses, function(x) mean(x[, "z.tl"], na.rm = T))
plot(sc, meanSES, pch = 21, bg = "forestgreen", xlab = "Уровень агрегации", ylab = "Средний SES MPD")
model <- lm(meanSES ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model, lwd = 2)
dev.off()

png("figures/report/herb-phylo-scaling-a.png", width = 600, height = 400)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
meanSES <- sapply(ha_sc_ses, function(x) mean(x[, "z.a.tl"], na.rm = T))
plot(sc, meanSES, pch = 21, bg = "forestgreen", xlab = "Уровень агрегации", ylab = "Средний SES MPD")
model <- lm(meanSES ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model, lwd = 2)
dev.off()

load("clean.data/source-pool-scaling.rda")

png("figures/report/herb-phylo-sp-scaling.png", width = 600, height = 400)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
fit <- lme(fixed = h.ses ~ r, random = ~ 1 | id, data = res, method = "REML", na.action = "na.omit")
ef <- effect("r", fit)
plot(h.ses ~ r, res, col = "forestgreen", pch = 16, xlab = "Размер видового пула", ylab = "SES MPD")
lines(ef$x$r, ef$fit, type = "l", lwd = 3, ylim = range(ef$lower, ef$upper), pch = 16)
lines(ef$x$r, ef$lower, lwd = 2, col = "grey40")
lines(ef$x$r, ef$upper, lwd = 2, col = "grey40")
dev.off()

png("figures/report/herb-phylo-sp-scaling-a.png", width = 600, height = 400)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
fit <- lme(fixed = h.ses.a ~ r, random = ~ 1 | id, data = res, method = "REML", na.action = "na.omit")
ef <- effect("r", fit)
plot(h.ses.a ~ r, res, col = "forestgreen", pch = 16, xlab = "Размер видового пула", ylab = "SES MPD")
lines(ef$x$r, ef$fit, type = "l", lwd = 3, ylim = range(ef$lower, ef$upper), pch = 16)
lines(ef$x$r, ef$lower, lwd = 2, col = "grey40")
lines(ef$x$r, ef$upper, lwd = 2, col = "grey40")
dev.off()

load("clean.data/coherent-scaling.rda")

png("figures/report/herb-phylo-coherent-scaling.png", width = 600, height = 400)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot(h.ses ~ sc, res, col = "forestgreen", pch = 16, xlab = "Уровень агрегации", ylab = "SES MPD")
if (anova(lm(h.ses ~ sc, res))[1, 5] < 0.05) abline(lm(h.ses ~ sc, res), lwd = 2)
dev.off()

png("figures/report/herb-phylo-coherent-scaling-a.png", width = 600, height = 400)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot(h.ses.a ~ sc, res, col = "forestgreen", pch = 16, xlab = "Уровень агрегации", ylab = "SES MPD")
if (anova(lm(h.ses.a ~ sc, res))[1, 5] < 0.05) abline(lm(h.ses.a ~ sc, res), lwd = 2)
dev.off()

