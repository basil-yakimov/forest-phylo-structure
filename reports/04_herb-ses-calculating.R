load('workspaces/herb-phylo.RData')

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
save.image("workspaces/herb-ses.RData")

#_____________________________________________________________________

load("workspaces/herb-ses.RData")

hgt <- read.table("data/hgt.txt")[[1]]

h.ses <- data.frame(hgt, taxa.labels = hses.tl$mpd.obs.z,
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

h.ses.p <- data.frame(hgt, taxa.labels = hses.tl$mpd.obs.p,
                      richness = hses.r$mpd.obs.p,
                      frequency = hses.f$mpd.obs.p,
                      sample.pool = hses.sp$mpd.obs.p,
                      phylogeny.pool = hses.pp$mpd.obs.p,
                      independentswap = hses.is$mpd.obs.p,
                      trialswap = hses.ts$mpd.obs.p,
                      taxa.labels.a = hses.a.tl$mpd.obs.p, 
                      richness.a = hses.a.r$mpd.obs.p,
                      frequency.a = hses.a.f$mpd.obs.p,
                      sample.pool.a = hses.a.sp$mpd.obs.p,
                      phylogeny.pool.a = hses.a.pp$mpd.obs.p,
                      independentswap.a = hses.a.is$mpd.obs.p,
                      trialswap.a = hses.a.ts$mpd.obs.p)
rm(list = ls(pattern = "hses"), dist.mat, tree)
save(h.ses, h.ses.p, file = "workspaces/herb-ses.rda")


