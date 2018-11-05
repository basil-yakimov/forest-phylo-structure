load('herb-phylo.RData')

library(ape)
library(picante)

dist.mat <- cophenetic(tree)

#_____________________________________________________________________

hses.tl <- ses.mpd(ha, dist.mat, null.model = "taxa.labels", runs = 999, iterations = 1000)

hses.r <- ses.mpd(ha, dist.mat, null.model = "richness", runs = 999, iterations = 1000)

hses.f <- ses.mpd(ha, dist.mat, null.model = "frequency", runs = 999, iterations = 1000)

hses.sp <- ses.mpd(ha, dist.mat, null.model = "sample.pool", runs = 999, iterations = 1000)

hses.pp <- ses.mpd(ha, dist.mat, null.model = "phylogeny.pool", runs = 999, iterations = 1000)

hses.is <- ses.mpd(ha, dist.mat, null.model = "independentswap", runs = 999, iterations = 1000)

hses.ts <- ses.mpd(ha, dist.mat, null.model = "trialswap", runs = 999, iterations = 1000)

#_____________________________________________________________________

hses.a.tl <- ses.mpd(ha, dist.mat, null.model = "taxa.labels", runs = 999, iterations = 1000, abundance.weighted = T)

hses.a.r <- ses.mpd(ha, dist.mat, null.model = "richness", runs = 999, iterations = 1000, abundance.weighted = T)

hses.a.f <- ses.mpd(ha, dist.mat, null.model = "frequency", runs = 999, iterations = 1000, abundance.weighted = T)

hses.a.sp <- ses.mpd(ha, dist.mat, null.model = "sample.pool", runs = 999, iterations = 1000, abundance.weighted = T)

hses.a.pp <- ses.mpd(ha, dist.mat, null.model = "phylogeny.pool", runs = 999, iterations = 1000, abundance.weighted = T)

hses.a.is <- ses.mpd(ha, dist.mat, null.model = "independentswap", runs = 999, iterations = 1000, abundance.weighted = T)

hses.a.ts <- ses.mpd(ha, dist.mat, null.model = "trialswap", runs = 999, iterations = 1000, abundance.weighted = T)

rm(my.add, my.bind.tip)
save.image("herb-phylo-analysis-ses.RData")

#_____________________________________________________________________

plot(1:96, hses.tl$mpd.obs)
plot(1:96, hses.tl$mpd.obs.z)

h.ses <- data.frame(n = 1:96, taxa.labels = hses.tl$mpd.obs.z,
                    richness = hses.r$mpd.obs.z,
                    frequency = hses.f$mpd.obs.z,
                    sample.pool = hses.sp$mpd.obs.z,
                    phylogeny.pool = hses.pp$mpd.obs.z,
                    independentswap = hses.is$mpd.obs.z,
                    trialswap = hses.ts$mpd.obs.z,
                    taxa.labels.a = hses.a.tl$mpd.obs.z, 
                    richness.a = hses.a.r$mpd.obs.z,
                    frequency.a = hses.a.f$mpd.obs.z,
                    sample.pool.a = hses.a.sp$mpd.obs.z,
                    phylogeny.pool.a = hses.a.pp$mpd.obs.z,
                    independentswap.a = hses.a.is$mpd.obs.z,
                    trialswap.a = hses.a.ts$mpd.obs.z)

hses.tl.m <- lm(taxa.labels ~ n, h.ses)
plot(h.ses$n, h.ses$taxa.labels)
abline(hses.tl.m)

shapiro.test(resid(hses.tl.m))
library(nortest)
lillie.test(resid(hses.tl.m))
ad.test(resid(hses.tl.m))

#_____________________________________________________________________

library(car)
ncvTest(hses.tl.m)

cor(h.ses$n, h.ses$taxa.labels)
cor.test(h.ses$n, h.ses$taxa.labels)

cor(h.ses$n, h.ses$taxa.labels, method = "spearman")
cor.test(h.ses$n, h.ses$taxa.labels, method = "spearman")

cor(h.ses$n, h.ses$taxa.labels, method = "kendall")
cor.test(h.ses$n, h.ses$taxa.labels, method = "kendall")

save.image("herb-phylo-analysis.RData")

