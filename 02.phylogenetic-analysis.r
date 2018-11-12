library(ape)
library(picante)

load("clean.data/wood-phylo.rda")
load("clean.data/herb-phylo.rda")

wood.mat <- cophenetic(wood.tree)
herb.mat <- cophenetic(herb.tree)

load("clean.data/tree-abund.rda")
load("clean.data/shrub-abund.rda")
load("clean.data/herb-abund.rda")

#---tree layer analysis ignoring abundances with various null models---#

t.ses.tl <- ses.mpd(ta, wood.mat, null.model = "taxa.labels", runs = 999, iterations = 1000)
t.ses.r  <- ses.mpd(ta, wood.mat, null.model = "richness", runs = 999, iterations = 1000)
t.ses.f  <- ses.mpd(ta, wood.mat, null.model = "frequency", runs = 999, iterations = 1000)
t.ses.sp <- ses.mpd(ta, wood.mat, null.model = "sample.pool", runs = 999, iterations = 1000)
t.ses.pp <- ses.mpd(ta, wood.mat, null.model = "phylogeny.pool", runs = 999, iterations = 1000)
t.ses.is <- ses.mpd(ta, wood.mat, null.model = "independentswap", runs = 999, iterations = 1000)
t.ses.ts <- ses.mpd(ta, wood.mat, null.model = "trialswap", runs = 999, iterations = 1000)

#---tree layer analysis considering abundances with various null models---#

t.ses.a.tl <- ses.mpd(ta, wood.mat, null.model = "taxa.labels", runs = 999, iterations = 1000, abundance.weighted = T)
t.ses.a.r  <- ses.mpd(ta, wood.mat, null.model = "richness", runs = 999, iterations = 1000, abundance.weighted = T)
t.ses.a.f  <- ses.mpd(ta, wood.mat, null.model = "frequency", runs = 999, iterations = 1000, abundance.weighted = T)
t.ses.a.sp <- ses.mpd(ta, wood.mat, null.model = "sample.pool", runs = 999, iterations = 1000, abundance.weighted = T)
t.ses.a.pp <- ses.mpd(ta, wood.mat, null.model = "phylogeny.pool", runs = 999, iterations = 1000, abundance.weighted = T)
t.ses.a.is <- ses.mpd(ta, wood.mat, null.model = "independentswap", runs = 999, iterations = 1000, abundance.weighted = T)
t.ses.a.ts <- ses.mpd(ta, wood.mat, null.model = "trialswap", runs = 999, iterations = 1000, abundance.weighted = T)

#---shrub layer analysis ignoring abundances with various null models---#

s.ses.tl <- ses.mpd(sa, wood.mat, null.model = "taxa.labels", runs = 999, iterations = 1000)
s.ses.r  <- ses.mpd(sa, wood.mat, null.model = "richness", runs = 999, iterations = 1000)
s.ses.f  <- ses.mpd(sa, wood.mat, null.model = "frequency", runs = 999, iterations = 1000)
s.ses.sp <- ses.mpd(sa, wood.mat, null.model = "sample.pool", runs = 999, iterations = 1000)
s.ses.pp <- ses.mpd(sa, wood.mat, null.model = "phylogeny.pool", runs = 999, iterations = 1000)
s.ses.is <- ses.mpd(sa, wood.mat, null.model = "independentswap", runs = 999, iterations = 1000)
s.ses.ts <- ses.mpd(sa, wood.mat, null.model = "trialswap", runs = 999, iterations = 1000)

#---tree layer analysis considering abundances with various null models---#

s.ses.a.tl <- ses.mpd(sa, wood.mat, null.model = "taxa.labels", runs = 999, iterations = 1000, abundance.weighted = T)
s.ses.a.r  <- ses.mpd(sa, wood.mat, null.model = "richness", runs = 999, iterations = 1000, abundance.weighted = T)
s.ses.a.f  <- ses.mpd(sa, wood.mat, null.model = "frequency", runs = 999, iterations = 1000, abundance.weighted = T)
s.ses.a.sp <- ses.mpd(sa, wood.mat, null.model = "sample.pool", runs = 999, iterations = 1000, abundance.weighted = T)
s.ses.a.pp <- ses.mpd(sa, wood.mat, null.model = "phylogeny.pool", runs = 999, iterations = 1000, abundance.weighted = T)
s.ses.a.is <- ses.mpd(sa, wood.mat, null.model = "independentswap", runs = 999, iterations = 1000, abundance.weighted = T)
s.ses.a.ts <- ses.mpd(sa, wood.mat, null.model = "trialswap", runs = 999, iterations = 1000, abundance.weighted = T)

#---herb layer analysis ignoring abundances with various null models---#

h.ses.tl <- ses.mpd(ha, herb.mat, null.model = "taxa.labels", runs = 999, iterations = 1000)
h.ses.r  <- ses.mpd(ha, herb.mat, null.model = "richness", runs = 999, iterations = 1000)
h.ses.f  <- ses.mpd(ha, herb.mat, null.model = "frequency", runs = 999, iterations = 1000)
h.ses.sp <- ses.mpd(ha, herb.mat, null.model = "sample.pool", runs = 999, iterations = 1000)
h.ses.pp <- ses.mpd(ha, herb.mat, null.model = "phylogeny.pool", runs = 999, iterations = 1000)
h.ses.is <- ses.mpd(ha, herb.mat, null.model = "independentswap", runs = 999, iterations = 1000)
h.ses.ts <- ses.mpd(ha, herb.mat, null.model = "trialswap", runs = 999, iterations = 1000)

#---herb layer analysis considering abundances with various null models---#

h.ses.a.tl <- ses.mpd(ha, herb.mat, null.model = "taxa.labels", runs = 999, iterations = 1000, abundance.weighted = T)
h.ses.a.r  <- ses.mpd(ha, herb.mat, null.model = "richness", runs = 999, iterations = 1000, abundance.weighted = T)
h.ses.a.f  <- ses.mpd(ha, herb.mat, null.model = "frequency", runs = 999, iterations = 1000, abundance.weighted = T)
h.ses.a.sp <- ses.mpd(ha, herb.mat, null.model = "sample.pool", runs = 999, iterations = 1000, abundance.weighted = T)
h.ses.a.pp <- ses.mpd(ha, herb.mat, null.model = "phylogeny.pool", runs = 999, iterations = 1000, abundance.weighted = T)
h.ses.a.is <- ses.mpd(ha, herb.mat, null.model = "independentswap", runs = 999, iterations = 1000, abundance.weighted = T)
h.ses.a.ts <- ses.mpd(ha, herb.mat, null.model = "trialswap", runs = 999, iterations = 1000, abundance.weighted = T)


#--- re-pack and save results ---#

hgt <- read.table("raw.data/hgt.txt")[[1]]

t.ses.z <- data.frame(hgt = hgt, 
                      taxa.labels = t.ses.tl$mpd.obs.z,
                      richness = t.ses.r$mpd.obs.z,
                      frequency = t.ses.f$mpd.obs.z,
                      sample.pool = t.ses.sp$mpd.obs.z,
                      phylogeny.pool = t.ses.pp$mpd.obs.z,
                      independentswap = t.ses.is$mpd.obs.z,
                      trialswap = t.ses.ts$mpd.obs.z,
                      taxa.labels.a = t.ses.a.tl$mpd.obs.z, 
                      richness.a = t.ses.a.r$mpd.obs.z,
                      frequency.a = t.ses.a.f$mpd.obs.z,
                      sample.pool.a = t.ses.a.sp$mpd.obs.z,
                      phylogeny.pool.a = t.ses.a.pp$mpd.obs.z,
                      independentswap.a = t.ses.a.is$mpd.obs.z,
                      trialswap.a = t.ses.a.ts$mpd.obs.z)

t.ses.p <- data.frame(hgt = hgt, 
                      taxa.labels = t.ses.tl$mpd.obs.p,
                      richness = t.ses.r$mpd.obs.p,
                      frequency = t.ses.f$mpd.obs.p,
                      sample.pool = t.ses.sp$mpd.obs.p,
                      phylogeny.pool = t.ses.pp$mpd.obs.p,
                      independentswap = t.ses.is$mpd.obs.p,
                      trialswap = t.ses.ts$mpd.obs.p,
                      taxa.labels.a = t.ses.a.tl$mpd.obs.p, 
                      richness.a = t.ses.a.r$mpd.obs.p,
                      frequency.a = t.ses.a.f$mpd.obs.p,
                      sample.pool.a = t.ses.a.sp$mpd.obs.p,
                      phylogeny.pool.a = t.ses.a.pp$mpd.obs.p,
                      independentswap.a = t.ses.a.is$mpd.obs.p,
                      trialswap.a = t.ses.a.ts$mpd.obs.p)

s.ses.z <- data.frame(hgt = hgt, 
                      taxa.labels = s.ses.tl$mpd.obs.z,
                      richness = s.ses.r$mpd.obs.z,
                      frequency = s.ses.f$mpd.obs.z,
                      sample.pool = s.ses.sp$mpd.obs.z,
                      phylogeny.pool = s.ses.pp$mpd.obs.z,
                      independentswap = s.ses.is$mpd.obs.z,
                      trialswap = s.ses.ts$mpd.obs.z,
                      taxa.labels.a = s.ses.a.tl$mpd.obs.z, 
                      richness.a = s.ses.a.r$mpd.obs.z,
                      frequency.a = s.ses.a.f$mpd.obs.z,
                      sample.pool.a = s.ses.a.sp$mpd.obs.z,
                      phylogeny.pool.a = s.ses.a.pp$mpd.obs.z,
                      independentswap.a = s.ses.a.is$mpd.obs.z,
                      trialswap.a = s.ses.a.ts$mpd.obs.z)

s.ses.p <- data.frame(hgt = hgt, 
                      taxa.labels = s.ses.tl$mpd.obs.p,
                      richness = s.ses.r$mpd.obs.p,
                      frequency = s.ses.f$mpd.obs.p,
                      sample.pool = s.ses.sp$mpd.obs.p,
                      phylogeny.pool = s.ses.pp$mpd.obs.p,
                      independentswap = s.ses.is$mpd.obs.p,
                      trialswap = s.ses.ts$mpd.obs.p,
                      taxa.labels.a = s.ses.a.tl$mpd.obs.p, 
                      richness.a = s.ses.a.r$mpd.obs.p,
                      frequency.a = s.ses.a.f$mpd.obs.p,
                      sample.pool.a = s.ses.a.sp$mpd.obs.p,
                      phylogeny.pool.a = s.ses.a.pp$mpd.obs.p,
                      independentswap.a = s.ses.a.is$mpd.obs.p,
                      trialswap.a = s.ses.a.ts$mpd.obs.p)

h.ses.z <- data.frame(hgt = hgt, 
                      taxa.labels = h.ses.tl$mpd.obs.z,
                      richness = h.ses.r$mpd.obs.z,
                      frequency = h.ses.f$mpd.obs.z,
                      sample.pool = h.ses.sp$mpd.obs.z,
                      phylogeny.pool = h.ses.pp$mpd.obs.z,
                      independentswap = h.ses.is$mpd.obs.z,
                      trialswap = h.ses.ts$mpd.obs.z,
                      taxa.labels.a = h.ses.a.tl$mpd.obs.z, 
                      richness.a = h.ses.a.r$mpd.obs.z,
                      frequency.a = h.ses.a.f$mpd.obs.z,
                      sample.pool.a = h.ses.a.sp$mpd.obs.z,
                      phylogeny.pool.a = h.ses.a.pp$mpd.obs.z,
                      independentswap.a = h.ses.a.is$mpd.obs.z,
                      trialswap.a = h.ses.a.ts$mpd.obs.z)

h.ses.p <- data.frame(hgt = hgt, 
                      taxa.labels = h.ses.tl$mpd.obs.p,
                      richness = h.ses.r$mpd.obs.p,
                      frequency = h.ses.f$mpd.obs.p,
                      sample.pool = h.ses.sp$mpd.obs.p,
                      phylogeny.pool = h.ses.pp$mpd.obs.p,
                      independentswap = h.ses.is$mpd.obs.p,
                      trialswap = h.ses.ts$mpd.obs.p,
                      taxa.labels.a = h.ses.a.tl$mpd.obs.p, 
                      richness.a = h.ses.a.r$mpd.obs.p,
                      frequency.a = h.ses.a.f$mpd.obs.p,
                      sample.pool.a = h.ses.a.sp$mpd.obs.p,
                      phylogeny.pool.a = h.ses.a.pp$mpd.obs.p,
                      independentswap.a = h.ses.a.is$mpd.obs.p,
                      trialswap.a = h.ses.a.ts$mpd.obs.p)

save(t.ses.z, t.ses.p, s.ses.z, s.ses.p, h.ses.z, h.ses.p, file = "clean.data/phylo-ses.rda")
