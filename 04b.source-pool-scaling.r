library(ape)
library(picante)
library(V.PhyloMaker)
library(dplyr)
library(readxl)


load("clean.data/tree-abund.rda")
load("clean.data/shrub-abund.rda")
load("clean.data/herb-abund.rda")

hgt <- read.table("raw.data/hgt.txt")[[1]]

load("clean.data/wood-phylo.rda")
load("clean.data/herb-phylo.rda")

ta <- ta[sample(nrow(ta)), ]
sa <- sa[sample(nrow(sa)), ]
ha <- ha[sample(nrow(ha)), ]

save(ta, file = "clean.data/tree-abund2.rda")
save(sa, file = "clean.data/shrub-abund2.rda")
save(ha, file = "clean.data/herb-abund2.rda")


ta <- ta[sample(nrow(ta)), ]
sa <- sa[sample(nrow(sa)), ]
ha <- ha[sample(nrow(ha)), ]

save(ta, file = "clean.data/tree-abund3.rda")
save(sa, file = "clean.data/shrub-abund3.rda")
save(ha, file = "clean.data/herb-abund3.rda")


#___________________________________________________________________


rand.mpd <- function(cdm, dist, ab = F)
{
  perm <- sample(nrow(dist))
  rownames(dist) <- colnames(dist) <- colnames(dist)[perm]
  mpd(cdm, dist, ab)
}


rmax <- floor((96-1)/2)
n <- sum(96 - 2*(1:rmax))

res <- data.frame(matrix(NA, nrow = n, ncol = 8))
colnames(res) <- c("r", "id", "t.ses", "t.ses.a", "s.ses", "s.ses.a", "h.ses", "h.ses.a")

t.list <- colnames(ta)
s.list <- colnames(sa)
h.list <- colnames(ha)

counter <- 1
for (r in 1:rmax)
{
  nloc <- 96 - 2*r
  for (jj in 1:nloc)
  {
    print(paste0("r = ", r, ", num = ", jj))
    
    res$r[counter]  <- r
    res$id[counter] <- r + jj
    
    # tree layer
    target <- ta[r + jj, ]
    present <- colSums(ta[jj : (jj + 2*r), ]) > 0
    target <- target[present]
    pool <- t.list[present]
    
    if (length(pool) > 1)
    {
      pool.dist <- cophenetic(drop.tip(wood.tree, wood.tree$tip.label[!(wood.tree$tip.label %in% pool)]))
      
      null.model <- replicate(500, rand.mpd(target, pool.dist))
      res$t.ses[counter] <- (mpd(target, pool.dist) - mean(null.model))/sd(null.model)
      
      null.model <- replicate(500, rand.mpd(target, pool.dist, T))
      res$t.ses.a[counter] <- (mpd(target, pool.dist, T) - mean(null.model))/sd(null.model)
    }
    
    # shrub layer
    target <- sa[r + jj, ]
    present <- colSums(sa[jj : (jj + 2*r), ]) > 0
    target <- target[present]
    pool <- s.list[present]
    
    if (length(pool) > 1)
    {
      pool.dist <- cophenetic(drop.tip(wood.tree, wood.tree$tip.label[!(wood.tree$tip.label %in% pool)]))
      
      null.model <- replicate(500, rand.mpd(target, pool.dist))
      res$s.ses[counter] <- (mpd(target, pool.dist) - mean(null.model))/sd(null.model)
      
      null.model <- replicate(500, rand.mpd(target, pool.dist, T))
      res$s.ses.a[counter] <- (mpd(target, pool.dist, T) - mean(null.model))/sd(null.model)
    }
    
    # herb layer
    target <- ha[r + jj, ]
    present <- colSums(ha[jj : (jj + 2*r), ]) > 0
    target <- target[present]
    pool <- h.list[present]
    
    if (length(pool) > 1)
    {
      pool.dist <- cophenetic(drop.tip(herb.tree, herb.tree$tip.label[!(herb.tree$tip.label %in% pool)]))
      
      null.model <- replicate(500, rand.mpd(target, pool.dist))
      res$h.ses[counter] <- (mpd(target, pool.dist) - mean(null.model))/sd(null.model)
      
      null.model <- replicate(500, rand.mpd(target, pool.dist, T))
      res$h.ses.a[counter] <- (mpd(target, pool.dist, T) - mean(null.model))/sd(null.model)
    }
    
    counter <- counter + 1
  }
}



save(res, file = "clean.data/source-pool-scaling.rda")
save(res, file = "clean.data/source-pool-scaling2.rda")
save(res, file = "clean.data/source-pool-scaling3.rda")

#_____________________________________________________________________

load("clean.data/source-pool-scaling.rda")

png("figures/source-pool-ses.png", 4000, 1500, pointsize = 75)
op <- par(mfrow = c(2,3), mar = c(3, 3, 2, 0.5))
plot(t.ses ~ r, res, col = "tomato", pch = 16)
plot(s.ses ~ r, res, col = "skyblue", pch = 16)
plot(h.ses ~ r, res, col = "forestgreen", pch = 16)

plot(1:rmax, tapply(res$t.ses, res$r, mean, na.rm = T), pch = 19, col = "tomato")
plot(1:rmax, tapply(res$s.ses, res$r, mean, na.rm = T), pch = 19, col = "skyblue")
plot(1:rmax, tapply(res$h.ses, res$r, mean, na.rm = T), pch = 19, col = "forestgreen")

dev.off()


png("figures/source-pool-ses-a.png", 4000, 1500, pointsize = 75)
op <- par(mfrow = c(2,3), mar = c(3, 3, 2, 0.5))
plot(t.ses.a ~ r, res, col = "tomato", pch = 16)
plot(s.ses.a ~ r, res, col = "skyblue", pch = 16)
plot(h.ses.a ~ r, res, col = "forestgreen", pch = 16)

plot(1:rmax, tapply(res$t.ses.a, res$r, mean, na.rm = T), pch = 19, col = "tomato")
plot(1:rmax, tapply(res$s.ses.a, res$r, mean, na.rm = T), pch = 19, col = "skyblue")
plot(1:rmax, tapply(res$h.ses.a, res$r, mean, na.rm = T), pch = 19, col = "forestgreen")

dev.off()

#_____________________________________________________________________

load("clean.data/source-pool-scaling2.rda")

png("figures/source-pool-ses2.png", 4000, 1500, pointsize = 75)
op <- par(mfrow = c(2,3), mar = c(3, 3, 2, 0.5))
plot(t.ses ~ r, res, col = "tomato", pch = 16)
plot(s.ses ~ r, res, col = "skyblue", pch = 16)
plot(h.ses ~ r, res, col = "forestgreen", pch = 16)

plot(1:rmax, tapply(res$t.ses, res$r, mean, na.rm = T), pch = 19, col = "tomato")
plot(1:rmax, tapply(res$s.ses, res$r, mean, na.rm = T), pch = 19, col = "skyblue")
plot(1:rmax, tapply(res$h.ses, res$r, mean, na.rm = T), pch = 19, col = "forestgreen")

dev.off()


png("figures/source-pool-ses-a2.png", 4000, 1500, pointsize = 75)
op <- par(mfrow = c(2,3), mar = c(3, 3, 2, 0.5))
plot(t.ses.a ~ r, res, col = "tomato", pch = 16)
plot(s.ses.a ~ r, res, col = "skyblue", pch = 16)
plot(h.ses.a ~ r, res, col = "forestgreen", pch = 16)

plot(1:rmax, tapply(res$t.ses.a, res$r, mean, na.rm = T), pch = 19, col = "tomato")
plot(1:rmax, tapply(res$s.ses.a, res$r, mean, na.rm = T), pch = 19, col = "skyblue")
plot(1:rmax, tapply(res$h.ses.a, res$r, mean, na.rm = T), pch = 19, col = "forestgreen")

dev.off()

#_____________________________________________________________________

load("clean.data/source-pool-scaling3.rda")

png("figures/source-pool-ses3.png", 4000, 1500, pointsize = 75)
op <- par(mfrow = c(2,3), mar = c(3, 3, 2, 0.5))
plot(t.ses ~ r, res, col = "tomato", pch = 16)
plot(s.ses ~ r, res, col = "skyblue", pch = 16)
plot(h.ses ~ r, res, col = "forestgreen", pch = 16)

plot(1:rmax, tapply(res$t.ses, res$r, mean, na.rm = T), pch = 19, col = "tomato")
plot(1:rmax, tapply(res$s.ses, res$r, mean, na.rm = T), pch = 19, col = "skyblue")
plot(1:rmax, tapply(res$h.ses, res$r, mean, na.rm = T), pch = 19, col = "forestgreen")

dev.off()


png("figures/source-pool-ses-a3.png", 4000, 1500, pointsize = 75)
op <- par(mfrow = c(2,3), mar = c(3, 3, 2, 0.5))
plot(t.ses.a ~ r, res, col = "tomato", pch = 16)
plot(s.ses.a ~ r, res, col = "skyblue", pch = 16)
plot(h.ses.a ~ r, res, col = "forestgreen", pch = 16)

plot(1:rmax, tapply(res$t.ses.a, res$r, mean, na.rm = T), pch = 19, col = "tomato")
plot(1:rmax, tapply(res$s.ses.a, res$r, mean, na.rm = T), pch = 19, col = "skyblue")
plot(1:rmax, tapply(res$h.ses.a, res$r, mean, na.rm = T), pch = 19, col = "forestgreen")

dev.off()



#load("clean.data/source-pool-scaling.rda")
#load("clean.data/source-pool-scaling2.rda")
#load("clean.data/source-pool-scaling3.rda")

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


op <- par(mfrow = c(3,2), mar = c(3, 3, 2, 0.5))

fit <- lme(fixed = t.ses.a ~ r, random = ~ 1 | id, data = res, method = "REML", na.action = "na.omit")
plot.ef(fit, col = "tomato")
title(main = coef(fit)[1, 2])

fit <- lm(t.ses.a ~ r + id, data = res)
plot.ef(fit, col = "tomato")
title(main = coef(fit)[2])


fit <- lme(fixed = s.ses.a ~ r, random = ~ 1 | id, data = res, method = "REML", na.action = "na.omit")
plot.ef(fit, col = "skyblue")
title(main = coef(fit)[1, 2])

fit <- lm(s.ses.a ~ r + id, data = res)
plot.ef(fit, col = "skyblue")
title(main = coef(fit)[2])

fit <- lme(fixed = h.ses.a ~ r, random = ~ 1 | id, data = res, method = "REML", na.action = "na.omit")
plot.ef(fit, col = "forestgreen")
title(main = coef(fit)[1, 2])

fit <- lm(h.ses.a ~ r + id, data = res)
plot.ef(fit, col = "forestgreen")
title(main = coef(fit)[2])

par(op)



op <- par(mfrow = c(3,2), mar = c(3, 3, 2, 0.5))

fit <- lme(fixed = t.ses ~ r, random = ~ 1 | id, data = res, method = "REML", na.action = "na.omit")
plot.ef(fit, col = "tomato")
title(main = coef(fit)[1, 2])

fit <- lm(t.ses ~ r + id, data = res)
plot.ef(fit, col = "tomato")
title(main = coef(fit)[2])


fit <- lme(fixed = s.ses ~ r, random = ~ 1 | id, data = res, method = "REML", na.action = "na.omit")
plot.ef(fit, col = "skyblue")
title(main = coef(fit)[1, 2])

fit <- lm(s.ses ~ r + id, data = res)
plot.ef(fit, col = "skyblue")
title(main = coef(fit)[2])

fit <- lme(fixed = h.ses ~ r, random = ~ 1 | id, data = res, method = "REML", na.action = "na.omit")
plot.ef(fit, col = "forestgreen")
title(main = coef(fit)[1, 2])

fit <- lm(h.ses ~ r + id, data = res)
plot.ef(fit, col = "forestgreen")
title(main = coef(fit)[2])

par(op)


op <- par(mfcol = c(3,2), mar = c(3, 3, 2, 0.5))

plot(0, 1, type = "n", xlim = c(0, 50), ylim = range(res$t.ses, na.rm = T))
for (ii in 2:99) lines(t.ses ~ r, data = res, subset = id == ii, col = "tomato")

plot(0, 1, type = "n", xlim = c(0, 50), ylim = range(res$s.ses, na.rm = T))
for (ii in 2:99) lines(s.ses ~ r, data = res, subset = id == ii, col = "skyblue")

plot(0, 1, type = "n", xlim = c(0, 50), ylim = range(res$h.ses, na.rm = T))
for (ii in 2:99) lines(h.ses ~ r, data = res, subset = id == ii, col = "forestgreen")

plot(0, 1, type = "n", xlim = c(0, 50), ylim = range(res$t.ses.a, na.rm = T))
for (ii in 2:99) lines(t.ses.a ~ r, data = res, subset = id == ii, col = "tomato")

plot(0, 1, type = "n", xlim = c(0, 50), ylim = range(res$s.ses.a, na.rm = T))
for (ii in 2:99) lines(s.ses.a ~ r, data = res, subset = id == ii, col = "skyblue")

plot(0, 1, type = "n", xlim = c(0, 50), ylim = range(res$h.ses.a, na.rm = T))
for (ii in 2:99) lines(h.ses.a ~ r, data = res, subset = id == ii, col = "forestgreen")

par(op)

