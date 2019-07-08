load("clean.data/beta.full.rda")

source("R/plot.ses.r")

png("figures/beta-tax-vs-dh.png", width = 3000, height = 2000)
op <- par(mfcol = c(3, 3), mar = c(4, 4, 0.5, 0.5), cex = 2)

plot.ses(beta$tt0, beta$dh, col = rgb(0.25, 0.75, 0.25, 0.3), lab = expression(beta[0]), xlab = expression(Delta * "hgt"))
plot.ses(beta$tt1, beta$dh, col = rgb(0.25, 0.75, 0.25, 0.3), lab = expression(beta[1]), xlab = expression(Delta * "hgt"))
plot.ses(beta$tt2, beta$dh, col = rgb(0.25, 0.75, 0.25, 0.3), lab = expression(beta[2]), xlab = expression(Delta * "hgt"))

plot.ses(beta$st0, beta$dh, col = rgb(0.25, 0.25, 0.75, 0.3), lab = expression(beta[0]), xlab = expression(Delta * "hgt"))
plot.ses(beta$st1, beta$dh, col = rgb(0.25, 0.25, 0.75, 0.3), lab = expression(beta[1]), xlab = expression(Delta * "hgt"))
plot.ses(beta$st2, beta$dh, col = rgb(0.25, 0.25, 0.75, 0.3), lab = expression(beta[2]), xlab = expression(Delta * "hgt"))

plot.ses(beta$ht0, beta$dh, col = rgb(0.75, 0.25, 0.25, 0.3), lab = expression(beta[0]), xlab = expression(Delta * "hgt"))
plot.ses(beta$ht1, beta$dh, col = rgb(0.75, 0.25, 0.25, 0.3), lab = expression(beta[1]), xlab = expression(Delta * "hgt"))
plot.ses(beta$ht2, beta$dh, col = rgb(0.75, 0.25, 0.25, 0.3), lab = expression(beta[2]), xlab = expression(Delta * "hgt"))

par(op)
dev.off()


png("figures/beta-tax-vs-mh.png", width = 3000, height = 2000)
op <- par(mfcol = c(3, 3), mar = c(4, 4, 0.5, 0.5), cex = 2)

plot.ses(beta$tt0, beta$mh, col = rgb(0.25, 0.75, 0.25, 0.3), lab = expression(beta[0]), xlab = "mean hgt")
plot.ses(beta$tt1, beta$mh, col = rgb(0.25, 0.75, 0.25, 0.3), lab = expression(beta[1]), xlab = "mean hgt")
plot.ses(beta$tt2, beta$mh, col = rgb(0.25, 0.75, 0.25, 0.3), lab = expression(beta[2]), xlab = "mean hgt")

plot.ses(beta$st0, beta$mh, col = rgb(0.25, 0.25, 0.75, 0.3), lab = expression(beta[0]), xlab = "mean hgt")
plot.ses(beta$st1, beta$mh, col = rgb(0.25, 0.25, 0.75, 0.3), lab = expression(beta[1]), xlab = "mean hgt")
plot.ses(beta$st2, beta$mh, col = rgb(0.25, 0.25, 0.75, 0.3), lab = expression(beta[2]), xlab = "mean hgt")

plot.ses(beta$ht0, beta$mh, col = rgb(0.75, 0.25, 0.25, 0.3), lab = expression(beta[0]), xlab = "mean hgt")
plot.ses(beta$ht1, beta$mh, col = rgb(0.75, 0.25, 0.25, 0.3), lab = expression(beta[1]), xlab = "mean hgt")
plot.ses(beta$ht2, beta$mh, col = rgb(0.75, 0.25, 0.25, 0.3), lab = expression(beta[2]), xlab = "mean hgt")

par(op)
dev.off()

#---#

png("figures/beta-phylo-vs-dh.png", width = 3000, height = 2000)
op <- par(mfcol = c(3, 3), mar = c(4, 4, 0.5, 0.5), cex = 2)

plot.ses(beta$tp0, beta$dh, col = rgb(0.25, 0.75, 0.25, 0.3), lab = expression("phylo - " * beta[0]), xlab = expression(Delta * "hgt"))
plot.ses(beta$tp1, beta$dh, col = rgb(0.25, 0.75, 0.25, 0.3), lab = expression("phylo - " * beta[1]), xlab = expression(Delta * "hgt"))
plot.ses(beta$tp2, beta$dh, col = rgb(0.25, 0.75, 0.25, 0.3), lab = expression("phylo - " * beta[2]), xlab = expression(Delta * "hgt"))

plot.ses(beta$sp0, beta$dh, col = rgb(0.25, 0.25, 0.75, 0.3), lab = expression("phylo - " * beta[0]), xlab = expression(Delta * "hgt"))
plot.ses(beta$sp1, beta$dh, col = rgb(0.25, 0.25, 0.75, 0.3), lab = expression("phylo - " * beta[1]), xlab = expression(Delta * "hgt"))
plot.ses(beta$sp2, beta$dh, col = rgb(0.25, 0.25, 0.75, 0.3), lab = expression("phylo - " * beta[2]), xlab = expression(Delta * "hgt"))

plot.ses(beta$hp0, beta$dh, col = rgb(0.75, 0.25, 0.25, 0.3), lab = expression("phylo - " * beta[0]), xlab = expression(Delta * "hgt"))
plot.ses(beta$hp1, beta$dh, col = rgb(0.75, 0.25, 0.25, 0.3), lab = expression("phylo - " * beta[1]), xlab = expression(Delta * "hgt"))
plot.ses(beta$hp2, beta$dh, col = rgb(0.75, 0.25, 0.25, 0.3), lab = expression("phylo - " * beta[2]), xlab = expression(Delta * "hgt"))

par(op)
dev.off()


png("figures/beta-phylo-vs-mh.png", width = 3000, height = 2000)
op <- par(mfcol = c(3, 3), mar = c(4, 4, 0.5, 0.5), cex = 2)

plot.ses(beta$tp0, beta$mh, col = rgb(0.25, 0.75, 0.25, 0.3), lab = expression("phylo - " * beta[0]), xlab = "mean hgt")
plot.ses(beta$tp1, beta$mh, col = rgb(0.25, 0.75, 0.25, 0.3), lab = expression("phylo - " * beta[1]), xlab = "mean hgt")
plot.ses(beta$tp2, beta$mh, col = rgb(0.25, 0.75, 0.25, 0.3), lab = expression("phylo - " * beta[2]), xlab = "mean hgt")

plot.ses(beta$sp0, beta$mh, col = rgb(0.25, 0.25, 0.75, 0.3), lab = expression("phylo - " * beta[0]), xlab = "mean hgt")
plot.ses(beta$sp1, beta$mh, col = rgb(0.25, 0.25, 0.75, 0.3), lab = expression("phylo - " * beta[1]), xlab = "mean hgt")
plot.ses(beta$sp2, beta$mh, col = rgb(0.25, 0.25, 0.75, 0.3), lab = expression("phylo - " * beta[2]), xlab = "mean hgt")

plot.ses(beta$hp0, beta$mh, col = rgb(0.75, 0.25, 0.25, 0.3), lab = expression("phylo - " * beta[0]), xlab = "mean hgt")
plot.ses(beta$hp1, beta$mh, col = rgb(0.75, 0.25, 0.25, 0.3), lab = expression("phylo - " * beta[1]), xlab = "mean hgt")
plot.ses(beta$hp2, beta$mh, col = rgb(0.75, 0.25, 0.25, 0.3), lab = expression("phylo - " * beta[2]), xlab = "mean hgt")

par(op)
dev.off()

#---#

png("figures/beta-phylo-vs-tax.png", width = 3000, height = 2000)
op <- par(mfcol = c(3, 3), mar = c(4, 4, 0.5, 0.5), cex = 2)

plot(beta$tp0 ~ beta$tt0, pch = 19, col = rgb(0.25, 0.75, 0.25, 0.3), ylab = expression("phylo - " * beta[0]), xlab = expression(beta[0]))
abline(0, 1)
points(beta$tp0.null ~ beta$tt0, pch = 19, col = rgb(0.01, 0.01, 0.01, 0.3))
fit <- lm(beta$tp0.null ~ beta$tt0)
abline(fit, lty = 2, lwd = 2, col = "grey")
beta$tp0.res <- beta$tp0 - predict(fit)

plot(beta$tp1 ~ beta$tt1, pch = 19, col = rgb(0.25, 0.75, 0.25, 0.3), ylab = expression("phylo - " * beta[1]), xlab = expression(beta[1]))
abline(0, 1)
points(beta$tp1.null ~ beta$tt1, pch = 19, col = rgb(0.01, 0.01, 0.01, 0.3))
fit <- lm(beta$tp1.null ~ beta$tt1)
abline(fit, lty = 2, lwd = 2, col = "grey")
beta$tp1.res <- beta$tp1 - predict(fit)

plot(beta$tp2 ~ beta$tt2, pch = 19, col = rgb(0.25, 0.75, 0.25, 0.3), ylab = expression("phylo - " * beta[2]), xlab = expression(beta[2]))
abline(0, 1)
points(beta$tp2.null ~ beta$tt2, pch = 19, col = rgb(0.01, 0.01, 0.01, 0.3))
fit <- lm(beta$tp2.null ~ beta$tt2)
abline(fit, lty = 2, lwd = 2, col = "grey")
beta$tp2.res <- beta$tp2 - predict(fit)

plot(beta$sp0 ~ beta$st0, pch = 19, col = rgb(0.25, 0.25, 0.75, 0.3), ylab = expression("phylo - " * beta[0]), xlab = expression(beta[0]))
abline(0, 1)
points(beta$sp0.null ~ beta$st0, pch = 19, col = rgb(0.01, 0.01, 0.01, 0.3))
fit <- lm(beta$sp0.null ~ beta$st0)
abline(fit, lty = 2, lwd = 2, col = "grey")
beta$sp0.res <- beta$sp0 - predict(fit)

plot(beta$sp1 ~ beta$st1, pch = 19, col = rgb(0.25, 0.25, 0.75, 0.3), ylab = expression("phylo - " * beta[1]), xlab = expression(beta[1]))
abline(0, 1)
points(beta$sp1.null ~ beta$st1, pch = 19, col = rgb(0.01, 0.01, 0.01, 0.3))
fit <- lm(beta$sp1.null ~ beta$st1)
abline(fit, lty = 2, lwd = 2, col = "grey")
beta$sp1.res <- beta$sp1 - predict(fit)

plot(beta$sp2 ~ beta$st2, pch = 19, col = rgb(0.25, 0.25, 0.75, 0.3), ylab = expression("phylo - " * beta[2]), xlab = expression(beta[2]))
abline(0, 1)
points(beta$sp2.null ~ beta$st2, pch = 19, col = rgb(0.01, 0.01, 0.01, 0.3))
fit <- lm(beta$sp2.null ~ beta$st2)
abline(fit, lty = 2, lwd = 2, col = "grey")
beta$sp2.res <- beta$sp2 - predict(fit)

plot(beta$hp0 ~ beta$ht0, pch = 19, col = rgb(0.75, 0.25, 0.25, 0.3), ylab = expression("phylo - " * beta[0]), xlab = expression(beta[0]))
abline(0, 1)
points(beta$hp0.null ~ beta$ht0, pch = 19, col = rgb(0.01, 0.01, 0.01, 0.3))
fit <- lm(beta$hp0.null ~ beta$ht0)
abline(fit, lty = 2, lwd = 2, col = "grey")
beta$hp0.res <- beta$hp0 - predict(fit)

plot(beta$hp1 ~ beta$ht1, pch = 19, col = rgb(0.75, 0.25, 0.25, 0.3), ylab = expression("phylo - " * beta[1]), xlab = expression(beta[1]))
abline(0, 1)
points(beta$hp1.null ~ beta$ht1, pch = 19, col = rgb(0.01, 0.01, 0.01, 0.3))
fit <- lm(beta$hp1.null ~ beta$ht1)
abline(fit, lty = 2, lwd = 2, col = "grey")
beta$hp1.res <- beta$hp1 - predict(fit)

plot(beta$hp2 ~ beta$ht2, pch = 19, col = rgb(0.75, 0.25, 0.25, 0.3), ylab = expression("phylo - " * beta[2]), xlab = expression(beta[2]))
abline(0, 1)
points(beta$hp2.null ~ beta$ht2, pch = 19, col = rgb(0.01, 0.01, 0.01, 0.3))
fit <- lm(beta$hp2.null ~ beta$ht2)
abline(fit, lty = 2, lwd = 2, col = "grey")
beta$hp2.res <- beta$hp2 - predict(fit)

par(op)
dev.off()


png("figures/beta-phylo-res-hist.png", width = 3000, height = 2000)
op <- par(mfcol = c(3, 3), mar = c(4, 4, 0.5, 0.5), cex = 2)

hist(beta$tp0.res, col = rgb(0.25, 0.75, 0.25, 0.3), xlab = expression("res - phylo - " * beta[0]), main = NULL)
abline(v = 0)
abline(v = mean(beta$tp0.res), lwd = 2, col = "grey")

hist(beta$tp1.res, col = rgb(0.25, 0.75, 0.25, 0.3), xlab = expression("res - phylo - " * beta[1]), main = NULL)
abline(v = 0)
abline(v = mean(beta$tp1.res), lwd = 2, col = "grey")

hist(beta$tp2.res, col = rgb(0.25, 0.75, 0.25, 0.3), xlab = expression("res - phylo - " * beta[2]), main = NULL)
abline(v = 0)
abline(v = mean(beta$tp2.res), lwd = 2, col = "grey")

hist(beta$sp0.res, col = rgb(0.25, 0.25, 0.75, 0.3), xlab = expression("res - phylo - " * beta[0]), main = NULL)
abline(v = 0)
abline(v = mean(beta$sp0.res), lwd = 2, col = "grey")

hist(beta$sp1.res, col = rgb(0.25, 0.25, 0.75, 0.3), xlab = expression("res - phylo - " * beta[1]), main = NULL)
abline(v = 0)
abline(v = mean(beta$sp1.res), lwd = 2, col = "grey")

hist(beta$sp2.res, col = rgb(0.25, 0.25, 0.75, 0.3), xlab = expression("res - phylo - " * beta[2]), main = NULL)
abline(v = 0)
abline(v = mean(beta$sp2.res), lwd = 2, col = "grey")

hist(beta$hp0.res, col = rgb(0.75, 0.25, 0.25, 0.3), xlab = expression("res - phylo - " * beta[0]), main = NULL)
abline(v = 0)
abline(v = mean(beta$hp0.res), lwd = 2, col = "grey")

hist(beta$hp1.res, col = rgb(0.75, 0.25, 0.25, 0.3), xlab = expression("res - phylo - " * beta[1]), main = NULL)
abline(v = 0)
abline(v = mean(beta$hp1.res), lwd = 2, col = "grey")

hist(beta$hp2.res, col = rgb(0.75, 0.25, 0.25, 0.3), xlab = expression("res - phylo - " * beta[2]), main = NULL)
abline(v = 0)
abline(v = mean(beta$hp2.res), lwd = 2, col = "grey")

par(op)
dev.off()


png("figures/beta-phylo-ses-hist.png", width = 3000, height = 2000)
op <- par(mfcol = c(3, 3), mar = c(4, 4, 0.5, 0.5), cex = 2)

hist(beta$tp0.ses, col = rgb(0.25, 0.75, 0.25, 0.3), xlab = expression("ses - phylo - " * beta[0]), main = NULL)
abline(v = 0)
abline(v = mean(beta$tp0.ses), lwd = 2, col = "grey")

hist(beta$tp1.ses, col = rgb(0.25, 0.75, 0.25, 0.3), xlab = expression("ses - phylo - " * beta[1]), main = NULL)
abline(v = 0)
abline(v = mean(beta$tp1.ses), lwd = 2, col = "grey")

hist(beta$tp2.ses, col = rgb(0.25, 0.75, 0.25, 0.3), xlab = expression("ses - phylo - " * beta[2]), main = NULL)
abline(v = 0)
abline(v = mean(beta$tp2.ses), lwd = 2, col = "grey")

hist(beta$sp0.ses, col = rgb(0.25, 0.25, 0.75, 0.3), xlab = expression("ses - phylo - " * beta[0]), main = NULL)
abline(v = 0)
abline(v = mean(beta$sp0.ses), lwd = 2, col = "grey")

hist(beta$sp1.ses, col = rgb(0.25, 0.25, 0.75, 0.3), xlab = expression("ses - phylo - " * beta[1]), main = NULL)
abline(v = 0)
abline(v = mean(beta$sp1.ses), lwd = 2, col = "grey")

hist(beta$sp2.ses, col = rgb(0.25, 0.25, 0.75, 0.3), xlab = expression("ses - phylo - " * beta[2]), main = NULL)
abline(v = 0)
abline(v = mean(beta$sp2.ses), lwd = 2, col = "grey")

hist(beta$hp0.ses, col = rgb(0.75, 0.25, 0.25, 0.3), xlab = expression("ses - phylo - " * beta[0]), main = NULL)
abline(v = 0)
abline(v = mean(beta$hp0.ses), lwd = 2, col = "grey")

hist(beta$hp1.ses, col = rgb(0.75, 0.25, 0.25, 0.3), xlab = expression("ses - phylo - " * beta[1]), main = NULL)
abline(v = 0)
abline(v = mean(beta$hp1.ses), lwd = 2, col = "grey")

hist(beta$hp2.ses, col = rgb(0.75, 0.25, 0.25, 0.3), xlab = expression("ses - phylo - " * beta[2]), main = NULL)
abline(v = 0)
abline(v = mean(beta$hp2.ses), lwd = 2, col = "grey")

par(op)
dev.off()


png("figures/beta-res-vs-ses.png", width = 3000, height = 2000)
op <- par(mfcol = c(3, 3), mar = c(4, 4, 0.5, 0.5), cex = 2)

plot.ses(beta$tp0.ses, beta$tp0.res, col = rgb(0.25, 0.75, 0.25, 0.3), lab = expression("res - phylo - " * beta[0]), xlab = expression("ses - phylo - " * beta[0]))
plot.ses(beta$tp1.ses, beta$tp1.res, col = rgb(0.25, 0.75, 0.25, 0.3), lab = expression("res - phylo - " * beta[1]), xlab = expression("ses - phylo - " * beta[1]))
plot.ses(beta$tp2.ses, beta$tp2.res, col = rgb(0.25, 0.75, 0.25, 0.3), lab = expression("res - phylo - " * beta[2]), xlab = expression("ses - phylo - " * beta[2]))

plot.ses(beta$sp0.ses, beta$sp0.res, col = rgb(0.25, 0.25, 0.75, 0.3), lab = expression("res - phylo - " * beta[0]), xlab = expression("ses - phylo - " * beta[0]))
plot.ses(beta$sp1.ses, beta$sp1.res, col = rgb(0.25, 0.25, 0.75, 0.3), lab = expression("res - phylo - " * beta[1]), xlab = expression("ses - phylo - " * beta[1]))
plot.ses(beta$sp2.ses, beta$sp2.res, col = rgb(0.25, 0.25, 0.75, 0.3), lab = expression("res - phylo - " * beta[2]), xlab = expression("ses - phylo - " * beta[2]))

plot.ses(beta$hp0.ses, beta$hp0.res, col = rgb(0.75, 0.25, 0.25, 0.3), lab = expression("res - phylo - " * beta[0]), xlab = expression("ses - phylo - " * beta[0]))
plot.ses(beta$hp1.ses, beta$hp1.res, col = rgb(0.75, 0.25, 0.25, 0.3), lab = expression("res - phylo - " * beta[1]), xlab = expression("ses - phylo - " * beta[1]))
plot.ses(beta$hp2.ses, beta$hp2.res, col = rgb(0.75, 0.25, 0.25, 0.3), lab = expression("res - phylo - " * beta[2]), xlab = expression("ses - phylo - " * beta[2]))

par(op)
dev.off()

#---#

png("figures/beta-res-vs-dh.png", width = 3000, height = 2000)
op <- par(mfcol = c(3, 3), mar = c(4, 4, 0.5, 0.5), cex = 2)

plot.ses(beta$tp0.res, beta$dh, col = rgb(0.25, 0.75, 0.25, 0.3), lab = expression("res - phylo - " * beta[0]), xlab = expression(Delta * "hgt"))
plot.ses(beta$tp1.res, beta$dh, col = rgb(0.25, 0.75, 0.25, 0.3), lab = expression("res - phylo - " * beta[1]), xlab = expression(Delta * "hgt"))
plot.ses(beta$tp2.res, beta$dh, col = rgb(0.25, 0.75, 0.25, 0.3), lab = expression("res - phylo - " * beta[2]), xlab = expression(Delta * "hgt"))

plot.ses(beta$sp0.res, beta$dh, col = rgb(0.25, 0.25, 0.75, 0.3), lab = expression("res - phylo - " * beta[0]), xlab = expression(Delta * "hgt"))
plot.ses(beta$sp1.res, beta$dh, col = rgb(0.25, 0.25, 0.75, 0.3), lab = expression("res - phylo - " * beta[1]), xlab = expression(Delta * "hgt"))
plot.ses(beta$sp2.res, beta$dh, col = rgb(0.25, 0.25, 0.75, 0.3), lab = expression("res - phylo - " * beta[2]), xlab = expression(Delta * "hgt"))

plot.ses(beta$hp0.res, beta$dh, col = rgb(0.75, 0.25, 0.25, 0.3), lab = expression("res - phylo - " * beta[0]), xlab = expression(Delta * "hgt"))
plot.ses(beta$hp1.res, beta$dh, col = rgb(0.75, 0.25, 0.25, 0.3), lab = expression("res - phylo - " * beta[1]), xlab = expression(Delta * "hgt"))
plot.ses(beta$hp2.res, beta$dh, col = rgb(0.75, 0.25, 0.25, 0.3), lab = expression("res - phylo - " * beta[2]), xlab = expression(Delta * "hgt"))

par(op)
dev.off()

png("figures/beta-res-vs-mh.png", width = 3000, height = 2000)
op <- par(mfcol = c(3, 3), mar = c(4, 4, 0.5, 0.5), cex = 2)

plot.ses(beta$tp0.res, beta$mh, col = rgb(0.25, 0.75, 0.25, 0.3), lab = expression("res - phylo - " * beta[0]), xlab = "mean hgt")
plot.ses(beta$tp1.res, beta$mh, col = rgb(0.25, 0.75, 0.25, 0.3), lab = expression("res - phylo - " * beta[1]), xlab = "mean hgt")
plot.ses(beta$tp2.res, beta$mh, col = rgb(0.25, 0.75, 0.25, 0.3), lab = expression("res - phylo - " * beta[2]), xlab = "mean hgt")

plot.ses(beta$sp0.res, beta$mh, col = rgb(0.25, 0.25, 0.75, 0.3), lab = expression("res - phylo - " * beta[0]), xlab = "mean hgt")
plot.ses(beta$sp1.res, beta$mh, col = rgb(0.25, 0.25, 0.75, 0.3), lab = expression("res - phylo - " * beta[1]), xlab = "mean hgt")
plot.ses(beta$sp2.res, beta$mh, col = rgb(0.25, 0.25, 0.75, 0.3), lab = expression("res - phylo - " * beta[2]), xlab = "mean hgt")

plot.ses(beta$hp0.res, beta$mh, col = rgb(0.75, 0.25, 0.25, 0.3), lab = expression("res - phylo - " * beta[0]), xlab = "mean hgt")
plot.ses(beta$hp1.res, beta$mh, col = rgb(0.75, 0.25, 0.25, 0.3), lab = expression("res - phylo - " * beta[1]), xlab = "mean hgt")
plot.ses(beta$hp2.res, beta$mh, col = rgb(0.75, 0.25, 0.25, 0.3), lab = expression("res - phylo - " * beta[2]), xlab = "mean hgt")

par(op)
dev.off()


png("figures/beta-ses-vs-dh.png", width = 3000, height = 2000)
op <- par(mfcol = c(3, 3), mar = c(4, 4, 0.5, 0.5), cex = 2)

plot.ses(beta$tp0.ses, beta$dh, col = rgb(0.25, 0.75, 0.25, 0.3), lab = expression("ses - phylo - " * beta[0]), xlab = expression(Delta * "hgt"))
plot.ses(beta$tp1.ses, beta$dh, col = rgb(0.25, 0.75, 0.25, 0.3), lab = expression("ses - phylo - " * beta[1]), xlab = expression(Delta * "hgt"))
plot.ses(beta$tp2.ses, beta$dh, col = rgb(0.25, 0.75, 0.25, 0.3), lab = expression("ses - phylo - " * beta[2]), xlab = expression(Delta * "hgt"))

plot.ses(beta$sp0.ses, beta$dh, col = rgb(0.25, 0.25, 0.75, 0.3), lab = expression("ses - phylo - " * beta[0]), xlab = expression(Delta * "hgt"))
plot.ses(beta$sp1.ses, beta$dh, col = rgb(0.25, 0.25, 0.75, 0.3), lab = expression("ses - phylo - " * beta[1]), xlab = expression(Delta * "hgt"))
plot.ses(beta$sp2.ses, beta$dh, col = rgb(0.25, 0.25, 0.75, 0.3), lab = expression("ses - phylo - " * beta[2]), xlab = expression(Delta * "hgt"))

plot.ses(beta$hp0.ses, beta$dh, col = rgb(0.75, 0.25, 0.25, 0.3), lab = expression("ses - phylo - " * beta[0]), xlab = expression(Delta * "hgt"))
plot.ses(beta$hp1.ses, beta$dh, col = rgb(0.75, 0.25, 0.25, 0.3), lab = expression("ses - phylo - " * beta[1]), xlab = expression(Delta * "hgt"))
plot.ses(beta$hp2.ses, beta$dh, col = rgb(0.75, 0.25, 0.25, 0.3), lab = expression("ses - phylo - " * beta[2]), xlab = expression(Delta * "hgt"))

par(op)
dev.off()


png("figures/beta-ses-vs-mh.png", width = 3000, height = 2000)
op <- par(mfcol = c(3, 3), mar = c(4, 4, 0.5, 0.5), cex = 2)

plot.ses(beta$tp0.ses, beta$mh, col = rgb(0.25, 0.75, 0.25, 0.3), lab = expression("ses - phylo - " * beta[0]), xlab = "mean hgt")
plot.ses(beta$tp1.ses, beta$mh, col = rgb(0.25, 0.75, 0.25, 0.3), lab = expression("ses - phylo - " * beta[1]), xlab = "mean hgt")
plot.ses(beta$tp2.ses, beta$mh, col = rgb(0.25, 0.75, 0.25, 0.3), lab = expression("ses - phylo - " * beta[2]), xlab = "mean hgt")

plot.ses(beta$sp0.ses, beta$mh, col = rgb(0.25, 0.25, 0.75, 0.3), lab = expression("ses - phylo - " * beta[0]), xlab = "mean hgt")
plot.ses(beta$sp1.ses, beta$mh, col = rgb(0.25, 0.25, 0.75, 0.3), lab = expression("ses - phylo - " * beta[1]), xlab = "mean hgt")
plot.ses(beta$sp2.ses, beta$mh, col = rgb(0.25, 0.25, 0.75, 0.3), lab = expression("ses - phylo - " * beta[2]), xlab = "mean hgt")

plot.ses(beta$hp0.ses, beta$mh, col = rgb(0.75, 0.25, 0.25, 0.3), lab = expression("ses - phylo - " * beta[0]), xlab = "mean hgt")
plot.ses(beta$hp1.ses, beta$mh, col = rgb(0.75, 0.25, 0.25, 0.3), lab = expression("ses - phylo - " * beta[1]), xlab = "mean hgt")
plot.ses(beta$hp2.ses, beta$mh, col = rgb(0.75, 0.25, 0.25, 0.3), lab = expression("ses - phylo - " * beta[2]), xlab = "mean hgt")

par(op)
dev.off()



#----------------------------------#

library(vegan)

my.as.dist <- function(x)
{
  dddd <- matrix(0, nrow = 96, ncol = 96)
  for (ii in 1:95) dddd[(ii+1):(96), ii] <- x[(sum(95:(95-ii+1))-(95-ii+1)+1):(sum(95:(95-ii+1)))]
  as.dist(dddd)
}


mantel(my.as.dist(beta$sp1), beta$st1)
mantel(my.as.dist(beta$sp1), beta$dh)

mantel.partial(my.as.dist(beta$sp1), beta$dh, beta$st1)

mantel(my.as.dist(beta$sp1.ses), beta$dh)

cor(beta$sp1.ses, beta$dh)



library(lmPerm)

fit <- lmp(hp2 ~ as.vector(ht2) + as.vector(dh) + as.vector(mh), beta)
summary(fit)

fit <- lmp(hp2 ~ as.vector(dh) * as.vector(mh), beta)
summary(fit)
