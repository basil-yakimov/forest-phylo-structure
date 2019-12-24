#library(ape)
library(picante)
#library(V.PhyloMaker)
#library(dplyr)
#library(readxl)


#load("clean.data/tree-abund.rda")
#load("clean.data/shrub-abund.rda")
#load("clean.data/herb-abund.rda")

#hgt <- read.table("raw.data/hgt.txt")[[1]]

load("clean.data/wood-phylo.rda")
load("clean.data/herb-phylo.rda")

#--------------------------------------------------------------#

load("clean.data/scaling.rda")

rand.mpd <- function(cdm, dist, ab = F)
{
  perm <- sample(nrow(dist))
  rownames(dist) <- colnames(dist) <- colnames(dist)[perm]
  mpd(cdm, dist, ab)
}

t.list <- colnames(ta_sc[[1]])
s.list <- colnames(sa_sc[[1]])
h.list <- colnames(ha_sc[[1]])


r <- 1

nums <- sapply(ta_sc, function(x) nrow(x))
n <- sum(nums - 2)

res <- data.frame(matrix(NA, nrow = n, ncol = 7))
colnames(res) <- c("sc", "t.ses", "t.ses.a", "s.ses", "s.ses.a", "h.ses", "h.ses.a")

counter <- 1

for (sc in 1:20)
{
  ta <- ta_sc[[sc]]
  sa <- sa_sc[[sc]]
  ha <- ha_sc[[sc]]
  
  nloc <- nrow(ta) - 2*1
  
  for (jj in 1:nloc)
  {
    print(paste0("sc = ", sc, ", num = ", jj))
    
    res$sc[counter]  <- sc
    #res$id[counter] <- r + jj
    
    # tree layer
    target <- ta[r + jj, , drop = F]
    present <- colSums(ta[jj : (jj + 2*r), ]) > 0
    target <- target[, present, drop = F]
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
    target <- sa[r + jj, , drop = F]
    present <- colSums(sa[jj : (jj + 2*r), ]) > 0
    target <- target[, present, drop = F]
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
    target <- ha[r + jj, , drop = F]
    present <- colSums(ha[jj : (jj + 2*r), ]) > 0
    target <- target[, present, drop = F]
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


save(res, file = "clean.data/coherent-scaling.rda")

#_____________________________________________________________________

load("clean.data/coherent-scaling.rda")
sc <- 1:20

png("figures/source-pool-ses.png", 4000, 1500, pointsize = 75)
op <- par(mfrow = c(2,3), mar = c(3, 3, 2, 0.5))
plot(t.ses ~ sc, res, col = "tomato", pch = 16)
abline(lm(t.ses ~ sc, res))
plot(s.ses ~ sc, res, col = "skyblue", pch = 16)
abline(lm(s.ses ~ sc, res))
plot(h.ses ~ sc, res, col = "forestgreen", pch = 16)
abline(lm(h.ses ~ sc, res))

plot(1:20, tapply(res$t.ses, res$sc, mean, na.rm = T), pch = 19, col = "tomato")
abline(lm(tapply(res$t.ses, res$sc, mean, na.rm = T) ~ sc))
plot(1:20, tapply(res$s.ses, res$sc, mean, na.rm = T), pch = 19, col = "skyblue")
abline(lm(tapply(res$s.ses, res$sc, mean, na.rm = T) ~ sc))
plot(1:20, tapply(res$h.ses, res$sc, mean, na.rm = T), pch = 19, col = "forestgreen")
abline(lm(tapply(res$h.ses, res$sc, mean, na.rm = T) ~ sc))
par(op)
dev.off()


png("figures/source-pool-ses-a.png", 4000, 1500, pointsize = 75)
op <- par(mfrow = c(2,3), mar = c(3, 3, 2, 0.5))
plot(t.ses.a ~ sc, res, col = "tomato", pch = 16)
abline(lm(t.ses.a ~ sc, res))
plot(s.ses.a ~ sc, res, col = "skyblue", pch = 16)
abline(lm(s.ses.a ~ sc, res))
plot(h.ses.a ~ sc, res, col = "forestgreen", pch = 16)
abline(lm(h.ses.a ~ sc, res))

plot(1:20, tapply(res$t.ses.a, res$sc, mean, na.rm = T), pch = 19, col = "tomato")
abline(lm(tapply(res$t.ses.a, res$sc, mean, na.rm = T) ~ sc))
plot(1:20, tapply(res$s.ses.a, res$sc, mean, na.rm = T), pch = 19, col = "skyblue")
abline(lm(tapply(res$s.ses.a, res$sc, mean, na.rm = T) ~ sc))
plot(1:20, tapply(res$h.ses.a, res$sc, mean, na.rm = T), pch = 19, col = "forestgreen")
abline(lm(tapply(res$h.ses.a, res$sc, mean, na.rm = T) ~ sc))
par(op)
dev.off()

