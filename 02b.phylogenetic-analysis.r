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

t.ses.tl <- ses.mntd(ta, wood.mat, null.model = "taxa.labels", runs = 999, iterations = 1000)
t.ses.r  <- ses.mntd(ta, wood.mat, null.model = "richness", runs = 999, iterations = 1000)
t.ses.f  <- ses.mntd(ta, wood.mat, null.model = "frequency", runs = 999, iterations = 1000)
t.ses.sp <- ses.mntd(ta, wood.mat, null.model = "sample.pool", runs = 999, iterations = 1000)
t.ses.pp <- ses.mntd(ta, wood.mat, null.model = "phylogeny.pool", runs = 999, iterations = 1000)
t.ses.is <- ses.mntd(ta, wood.mat, null.model = "independentswap", runs = 999, iterations = 1000)
t.ses.ts <- ses.mntd(ta, wood.mat, null.model = "trialswap", runs = 999, iterations = 1000)

#---tree layer analysis considering abundances with various null models---#

t.ses.a.tl <- ses.mntd(ta, wood.mat, null.model = "taxa.labels", runs = 999, iterations = 1000, abundance.weighted = T)
t.ses.a.r  <- ses.mntd(ta, wood.mat, null.model = "richness", runs = 999, iterations = 1000, abundance.weighted = T)
t.ses.a.f  <- ses.mntd(ta, wood.mat, null.model = "frequency", runs = 999, iterations = 1000, abundance.weighted = T)
t.ses.a.sp <- ses.mntd(ta, wood.mat, null.model = "sample.pool", runs = 999, iterations = 1000, abundance.weighted = T)
t.ses.a.pp <- ses.mntd(ta, wood.mat, null.model = "phylogeny.pool", runs = 999, iterations = 1000, abundance.weighted = T)
t.ses.a.is <- ses.mntd(ta, wood.mat, null.model = "independentswap", runs = 999, iterations = 1000, abundance.weighted = T)
t.ses.a.ts <- ses.mntd(ta, wood.mat, null.model = "trialswap", runs = 999, iterations = 1000, abundance.weighted = T)

#---shrub layer analysis ignoring abundances with various null models---#

s.ses.tl <- ses.mntd(sa, wood.mat, null.model = "taxa.labels", runs = 999, iterations = 1000)
s.ses.r  <- ses.mntd(sa, wood.mat, null.model = "richness", runs = 999, iterations = 1000)
s.ses.f  <- ses.mntd(sa, wood.mat, null.model = "frequency", runs = 999, iterations = 1000)
s.ses.sp <- ses.mntd(sa, wood.mat, null.model = "sample.pool", runs = 999, iterations = 1000)
s.ses.pp <- ses.mntd(sa, wood.mat, null.model = "phylogeny.pool", runs = 999, iterations = 1000)
s.ses.is <- ses.mntd(sa, wood.mat, null.model = "independentswap", runs = 999, iterations = 1000)
s.ses.ts <- ses.mntd(sa, wood.mat, null.model = "trialswap", runs = 999, iterations = 1000)

#---tree layer analysis considering abundances with various null models---#

s.ses.a.tl <- ses.mntd(sa, wood.mat, null.model = "taxa.labels", runs = 999, iterations = 1000, abundance.weighted = T)
s.ses.a.r  <- ses.mntd(sa, wood.mat, null.model = "richness", runs = 999, iterations = 1000, abundance.weighted = T)
s.ses.a.f  <- ses.mntd(sa, wood.mat, null.model = "frequency", runs = 999, iterations = 1000, abundance.weighted = T)
s.ses.a.sp <- ses.mntd(sa, wood.mat, null.model = "sample.pool", runs = 999, iterations = 1000, abundance.weighted = T)
s.ses.a.pp <- ses.mntd(sa, wood.mat, null.model = "phylogeny.pool", runs = 999, iterations = 1000, abundance.weighted = T)
s.ses.a.is <- ses.mntd(sa, wood.mat, null.model = "independentswap", runs = 999, iterations = 1000, abundance.weighted = T)
s.ses.a.ts <- ses.mntd(sa, wood.mat, null.model = "trialswap", runs = 999, iterations = 1000, abundance.weighted = T)

#---herb layer analysis ignoring abundances with various null models---#

h.ses.tl <- ses.mntd(ha, herb.mat, null.model = "taxa.labels", runs = 999, iterations = 1000)
h.ses.r  <- ses.mntd(ha, herb.mat, null.model = "richness", runs = 999, iterations = 1000)
h.ses.f  <- ses.mntd(ha, herb.mat, null.model = "frequency", runs = 999, iterations = 1000)
h.ses.sp <- ses.mntd(ha, herb.mat, null.model = "sample.pool", runs = 999, iterations = 1000)
h.ses.pp <- ses.mntd(ha, herb.mat, null.model = "phylogeny.pool", runs = 999, iterations = 1000)
h.ses.is <- ses.mntd(ha, herb.mat, null.model = "independentswap", runs = 999, iterations = 1000)
h.ses.ts <- ses.mntd(ha, herb.mat, null.model = "trialswap", runs = 999, iterations = 1000)

#---herb layer analysis considering abundances with various null models---#

h.ses.a.tl <- ses.mntd(ha, herb.mat, null.model = "taxa.labels", runs = 999, iterations = 1000, abundance.weighted = T)
h.ses.a.r  <- ses.mntd(ha, herb.mat, null.model = "richness", runs = 999, iterations = 1000, abundance.weighted = T)
h.ses.a.f  <- ses.mntd(ha, herb.mat, null.model = "frequency", runs = 999, iterations = 1000, abundance.weighted = T)
h.ses.a.sp <- ses.mntd(ha, herb.mat, null.model = "sample.pool", runs = 999, iterations = 1000, abundance.weighted = T)
h.ses.a.pp <- ses.mntd(ha, herb.mat, null.model = "phylogeny.pool", runs = 999, iterations = 1000, abundance.weighted = T)
h.ses.a.is <- ses.mntd(ha, herb.mat, null.model = "independentswap", runs = 999, iterations = 1000, abundance.weighted = T)
h.ses.a.ts <- ses.mntd(ha, herb.mat, null.model = "trialswap", runs = 999, iterations = 1000, abundance.weighted = T)


#--- re-pack and save results ---#

hgt <- read.table("raw.data/hgt.txt")[[1]]

t.ses.z <- data.frame(hgt = hgt, 
                      taxa.labels = t.ses.tl$mntd.obs.z,
                      richness = t.ses.r$mntd.obs.z,
                      frequency = t.ses.f$mntd.obs.z,
                      sample.pool = t.ses.sp$mntd.obs.z,
                      phylogeny.pool = t.ses.pp$mntd.obs.z,
                      independentswap = t.ses.is$mntd.obs.z,
                      trialswap = t.ses.ts$mntd.obs.z,
                      taxa.labels.a = t.ses.a.tl$mntd.obs.z, 
                      richness.a = t.ses.a.r$mntd.obs.z,
                      frequency.a = t.ses.a.f$mntd.obs.z,
                      sample.pool.a = t.ses.a.sp$mntd.obs.z,
                      phylogeny.pool.a = t.ses.a.pp$mntd.obs.z,
                      independentswap.a = t.ses.a.is$mntd.obs.z,
                      trialswap.a = t.ses.a.ts$mntd.obs.z)

t.ses.p <- data.frame(hgt = hgt, 
                      taxa.labels = t.ses.tl$mntd.obs.p,
                      richness = t.ses.r$mntd.obs.p,
                      frequency = t.ses.f$mntd.obs.p,
                      sample.pool = t.ses.sp$mntd.obs.p,
                      phylogeny.pool = t.ses.pp$mntd.obs.p,
                      independentswap = t.ses.is$mntd.obs.p,
                      trialswap = t.ses.ts$mntd.obs.p,
                      taxa.labels.a = t.ses.a.tl$mntd.obs.p, 
                      richness.a = t.ses.a.r$mntd.obs.p,
                      frequency.a = t.ses.a.f$mntd.obs.p,
                      sample.pool.a = t.ses.a.sp$mntd.obs.p,
                      phylogeny.pool.a = t.ses.a.pp$mntd.obs.p,
                      independentswap.a = t.ses.a.is$mntd.obs.p,
                      trialswap.a = t.ses.a.ts$mntd.obs.p)

s.ses.z <- data.frame(hgt = hgt, 
                      taxa.labels = s.ses.tl$mntd.obs.z,
                      richness = s.ses.r$mntd.obs.z,
                      frequency = s.ses.f$mntd.obs.z,
                      sample.pool = s.ses.sp$mntd.obs.z,
                      phylogeny.pool = s.ses.pp$mntd.obs.z,
                      independentswap = s.ses.is$mntd.obs.z,
                      trialswap = s.ses.ts$mntd.obs.z,
                      taxa.labels.a = s.ses.a.tl$mntd.obs.z, 
                      richness.a = s.ses.a.r$mntd.obs.z,
                      frequency.a = s.ses.a.f$mntd.obs.z,
                      sample.pool.a = s.ses.a.sp$mntd.obs.z,
                      phylogeny.pool.a = s.ses.a.pp$mntd.obs.z,
                      independentswap.a = s.ses.a.is$mntd.obs.z,
                      trialswap.a = s.ses.a.ts$mntd.obs.z)

s.ses.p <- data.frame(hgt = hgt, 
                      taxa.labels = s.ses.tl$mntd.obs.p,
                      richness = s.ses.r$mntd.obs.p,
                      frequency = s.ses.f$mntd.obs.p,
                      sample.pool = s.ses.sp$mntd.obs.p,
                      phylogeny.pool = s.ses.pp$mntd.obs.p,
                      independentswap = s.ses.is$mntd.obs.p,
                      trialswap = s.ses.ts$mntd.obs.p,
                      taxa.labels.a = s.ses.a.tl$mntd.obs.p, 
                      richness.a = s.ses.a.r$mntd.obs.p,
                      frequency.a = s.ses.a.f$mntd.obs.p,
                      sample.pool.a = s.ses.a.sp$mntd.obs.p,
                      phylogeny.pool.a = s.ses.a.pp$mntd.obs.p,
                      independentswap.a = s.ses.a.is$mntd.obs.p,
                      trialswap.a = s.ses.a.ts$mntd.obs.p)

h.ses.z <- data.frame(hgt = hgt, 
                      taxa.labels = h.ses.tl$mntd.obs.z,
                      richness = h.ses.r$mntd.obs.z,
                      frequency = h.ses.f$mntd.obs.z,
                      sample.pool = h.ses.sp$mntd.obs.z,
                      phylogeny.pool = h.ses.pp$mntd.obs.z,
                      independentswap = h.ses.is$mntd.obs.z,
                      trialswap = h.ses.ts$mntd.obs.z,
                      taxa.labels.a = h.ses.a.tl$mntd.obs.z, 
                      richness.a = h.ses.a.r$mntd.obs.z,
                      frequency.a = h.ses.a.f$mntd.obs.z,
                      sample.pool.a = h.ses.a.sp$mntd.obs.z,
                      phylogeny.pool.a = h.ses.a.pp$mntd.obs.z,
                      independentswap.a = h.ses.a.is$mntd.obs.z,
                      trialswap.a = h.ses.a.ts$mntd.obs.z)

h.ses.p <- data.frame(hgt = hgt, 
                      taxa.labels = h.ses.tl$mntd.obs.p,
                      richness = h.ses.r$mntd.obs.p,
                      frequency = h.ses.f$mntd.obs.p,
                      sample.pool = h.ses.sp$mntd.obs.p,
                      phylogeny.pool = h.ses.pp$mntd.obs.p,
                      independentswap = h.ses.is$mntd.obs.p,
                      trialswap = h.ses.ts$mntd.obs.p,
                      taxa.labels.a = h.ses.a.tl$mntd.obs.p, 
                      richness.a = h.ses.a.r$mntd.obs.p,
                      frequency.a = h.ses.a.f$mntd.obs.p,
                      sample.pool.a = h.ses.a.sp$mntd.obs.p,
                      phylogeny.pool.a = h.ses.a.pp$mntd.obs.p,
                      independentswap.a = h.ses.a.is$mntd.obs.p,
                      trialswap.a = h.ses.a.ts$mntd.obs.p)

save(t.ses.z, t.ses.p, s.ses.z, s.ses.p, h.ses.z, h.ses.p, file = "clean.data/phylo-ses-mntd.rda")
