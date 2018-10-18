load('tree-shrub-phylo.RData')

library(ape)
dist.mat <- cophenetic(tree)              
library(picante)
mpd.mat <- mpd(ta, cophenetic(tree), abundance.weighted = F) #NA
smpd.mat <- mpd(sa, cophenetic(tree), abundance.weighted = F)

#_____________________________________________________________________

mpd.mat.a <- mpd(ta, cophenetic(tree), abundance.weighted = T)
smpd.mat.a <- mpd(sa, cophenetic(tree), abundance.weighted = T)

#_____________________________________________________________________

ses <- ses.mpd(ta, dist.mat, null.model = "taxa.labels", runs = 999, iterations = 1000)

plot(1:96, ses$mpd.obs)
plot(1:96, ses$mpd.obs.z)
plot(1:96, ses$mpd.obs.p)

ses.a <- ses.mpd(sa, dist.mat, null.model = "taxa.labels", runs = 999, iterations = 1000)
plot(1:96, ses$mpd.obs)
plot(1:96, ses$mpd.obs.z)
plot(1:96, ses$mpd.obs.p)

