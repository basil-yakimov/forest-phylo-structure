load('tree-shrub-phylo.RData')

library(ape)
dist.mat <- cophenetic(tree)              
library(picante)
mpd.mat <- mpd(ta, dist.mat, abundance.weighted = F)
smpd.mat <- mpd(sa, dist.mat, abundance.weighted = F)

#_____________________________________________________________________

mpd.mat.a <- mpd(ta, dist.mat, abundance.weighted = T)
smpd.mat.a <- mpd(sa, dist.mat, abundance.weighted = T)

#_____________________________________________________________________

tses.tl <- ses.mpd(ta, dist.mat, null.model = "taxa.labels", runs = 999, iterations = 1000)

tses.r <- ses.mpd(ta, dist.mat, null.model = "richness", runs = 999, iterations = 1000)

tses.f <- ses.mpd(ta, dist.mat, null.model = "frequency", runs = 999, iterations = 1000)

tses.sp <- ses.mpd(ta, dist.mat, null.model = "sample.pool", runs = 999, iterations = 1000)

tses.pp <- ses.mpd(ta, dist.mat, null.model = "phylogeny.pool", runs = 999, iterations = 1000)

tses.is <- ses.mpd(ta, dist.mat, null.model = "independentswap", runs = 999, iterations = 1000)

tses.ts <- ses.mpd(ta, dist.mat, null.model = "trialswap", runs = 999, iterations = 1000)

#_____________________________________________________________________

tses.a.tl <- ses.mpd(ta, dist.mat, null.model = "taxa.labels", runs = 999, iterations = 1000, abundance.weighted = T)

tses.a.r <- ses.mpd(ta, dist.mat, null.model = "richness", runs = 999, iterations = 1000, abundance.weighted = T)

tses.a.f <- ses.mpd(ta, dist.mat, null.model = "frequency", runs = 999, iterations = 1000, abundance.weighted = T)

tses.a.sp <- ses.mpd(ta, dist.mat, null.model = "sample.pool", runs = 999, iterations = 1000, abundance.weighted = T)

tses.a.pp <- ses.mpd(ta, dist.mat, null.model = "phylogeny.pool", runs = 999, iterations = 1000, abundance.weighted = T)

tses.a.is <- ses.mpd(ta, dist.mat, null.model = "independentswap", runs = 999, iterations = 1000, abundance.weighted = T)

tses.a.ts  <- ses.mpd(ta, dist.mat, null.model = "trialswap", runs = 999, iterations = 1000, abundance.weighted = T)

#_____________________________________________________________________

sses.tl <- ses.mpd(sa, dist.mat, null.model = "taxa.labels", runs = 999, iterations = 1000)

sses.r <- ses.mpd(sa, dist.mat, null.model = "richness", runs = 999, iterations = 1000)

sses.f <- ses.mpd(sa, dist.mat, null.model = "frequency", runs = 999, iterations = 1000)

sses.sp <- ses.mpd(sa, dist.mat, null.model = "sample.pool", runs = 999, iterations = 1000)

sses.pp <- ses.mpd(sa, dist.mat, null.model = "phylogeny.pool", runs = 999, iterations = 1000)

sses.is <- ses.mpd(sa, dist.mat, null.model = "independentswap", runs = 999, iterations = 1000)

sses.ts <- ses.mpd(sa, dist.mat, null.model = "trialswap", runs = 999, iterations = 1000)

#_____________________________________________________________________

sses.a.tl <- ses.mpd(sa, dist.mat, null.model = "taxa.labels", runs = 999, iterations = 1000, abundance.weighted = T)

sses.a.r <- ses.mpd(sa, dist.mat, null.model = "richness", runs = 999, iterations = 1000, abundance.weighted = T)

sses.a.f <- ses.mpd(sa, dist.mat, null.model = "frequency", runs = 999, iterations = 1000, abundance.weighted = T)

sses.a.sp <- ses.mpd(sa, dist.mat, null.model = "sample.pool", runs = 999, iterations = 1000, abundance.weighted = T)

sses.a.pp <- ses.mpd(sa, dist.mat, null.model = "phylogeny.pool", runs = 999, iterations = 1000, abundance.weighted = T)

sses.a.is <- ses.mpd(sa, dist.mat, null.model = "independentswap", runs = 999, iterations = 1000, abundance.weighted = T)

sses.a.ts <- ses.mpd(sa, dist.mat, null.model = "trialswap", runs = 999, iterations = 1000, abundance.weighted = T)

rm(mpd.mat, mpd.mat.a, smpd.mat, smpd.mat.a, my.add, my.bind.tip)
save.image("tree-shrub-phylo-analysis-ses.RData")

#_____________________________________________________________________

plot(1:96, tses.tl$mpd.obs)
plot(1:96, tses.tl$mpd.obs.z)


hgt <- read.table("hgt.txt")[[1]]
  
  
  
t.ses <- data.frame(hgt, taxa.labels = tses.tl$mpd.obs.z,
                richness = tses.r$mpd.obs.z,
                frequency = tses.f$mpd.obs.z,
                sample.pool = tses.sp$mpd.obs.z,
                phylogeny.pool = tses.pp$mpd.obs.z,
                independentswap = tses.is$mpd.obs.z,
                trialswap = tses.ts$mpd.obs.z,
                taxa.labels.a = tses.a.tl$mpd.obs.z, 
                richness.a = tses.a.r$mpd.obs.z,
                frequency.a = tses.a.f$mpd.obs.z,
                sample.pool.a = tses.a.sp$mpd.obs.z,
                phylogeny.pool.a = tses.a.pp$mpd.obs.z,
                independentswap.a = tses.a.is$mpd.obs.z,
                trialswap.a = tses.a.ts$mpd.obs.z)

s.ses <- data.frame(hgt, taxa.labels = sses.tl$mpd.obs.z,
                    richness = sses.r$mpd.obs.z,
                    frequency = sses.f$mpd.obs.z,
                    sample.pool = sses.sp$mpd.obs.z,
                    phylogeny.pool = sses.pp$mpd.obs.z,
                    independentswap = sses.is$mpd.obs.z,
                    trialswap = sses.ts$mpd.obs.z,
                    taxa.labels.a = sses.a.tl$mpd.obs.z, 
                    richness.a = sses.a.r$mpd.obs.z,
                    frequency.a = sses.a.f$mpd.obs.z,
                    sample.pool.a = sses.a.sp$mpd.obs.z,
                    phylogeny.pool.a = sses.a.pp$mpd.obs.z,
                    independentswap.a = sses.a.is$mpd.obs.z,
                    trialswap.a = sses.a.ts$mpd.obs.z)

tses.tl.m <- lm(taxa.labels ~ hgt, t.ses)
plot(t.ses$hgt, t.ses$taxa.labels, pch = 21, bg = "steelblue")
abline(tses.tl.m)
title(main = paste0("taxa.labels - ", round(sp$estimate, dig = 3), " - ", round(sp$p.value, 3)))



shapiro.test(resid(tses.tl.m))
library(nortest)
lillie.test(resid(tses.tl.m))
ad.test(resid(tses.tl.m))

library(car)
ncvTest(tses.tl.m)
summary(tses.tl.m)

cor.test(t.ses$hgt[!is.na(t.ses$taxa.labels)], t.ses$taxa.labels[!is.na(t.ses$taxa.labels)], method = "spearman")


#_____________________________________________________________________

cor(t.ses$n[!is.na(t.ses$taxa.labels)], t.ses$taxa.labels[!is.na(t.ses$taxa.labels)])
cor.test(t.ses$n[!is.na(t.ses$taxa.labels)], t.ses$taxa.labels[!is.na(t.ses$taxa.labels)])

cor(t.ses$n[!is.na(t.ses$taxa.labels)], t.ses$taxa.labels[!is.na(t.ses$taxa.labels)], method = "spearman")
sp <- cor.test(t.ses$hgt[!is.na(t.ses$taxa.labels)], t.ses$taxa.labels[!is.na(t.ses$taxa.labels)], method = "spearman")

cor(t.ses$n[!is.na(t.ses$taxa.labels)], t.ses$taxa.labels[!is.na(t.ses$taxa.labels)], method = "kendall")
cor.test(t.ses$n[!is.na(t.ses$taxa.labels)], t.ses$taxa.labels[!is.na(t.ses$taxa.labels)], method = "kendall")

#_____________________________________________________________________

cor(s.ses$n, s.ses$taxa.labels)
cor.test(s.ses$n, s.ses$taxa.labels)

cor(s.ses$n, s.ses$taxa.labels, method = "spearman")
cor.test(s.ses$n, s.ses$taxa.labels, method = "spearman")

cor(s.ses$n, s.ses$taxa.labels, method = "kendall")
cor.test(s.ses$n, s.ses$taxa.labels, method = "kendall")

save.image("tree-shrub-phylo-analysis.RData")
