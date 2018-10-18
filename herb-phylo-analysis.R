load('herb-phylo.RData')

library(ape)
library(picante)
mpd.mat <- mpd(ha, cophenetic(tree), abundance.weighted = F)

#_____________________________________________________________________

mpd.mat.a <- mpd(ha, cophenetic(tree), abundance.weighted = T)

#_____________________________________________________________________

ses <- ses.mpd(ha, dist.mat, null.model = "taxa.labels", runs = 999, iterations = 1000)

plot(1:96, ses$mpd.obs)
plot(1:96, ses$mpd.obs.z)
plot(1:96, ses$mpd.obs.p)
