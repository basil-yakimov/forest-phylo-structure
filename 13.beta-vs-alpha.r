library(vegan)
library(picante)

source("R/plot.ses.r")

load("clean.data/wood-phylo.rda")
load("clean.data/herb-phylo.rda")

wood.mat <- cophenetic(wood.tree)
herb.mat <- cophenetic(herb.tree)

load("clean.data/tree-abund.rda")
load("clean.data/shrub-abund.rda")
load("clean.data/herb-abund.rda")

hgt <- read.table("raw.data/hgt.txt")[[1]]


#---#

load("clean.data/phylo-ses-pqd.rda")

#---#

load("clean.data/pw-beta-tl.rda")


op <- par(mfcol = c(2, 3))
plot.ses(beta.ses$td0[c(T, F)], (ta.ses$mpd[seq(1, 95, by = 2)] + ta.ses$mpd[seq(2, 96, by = 2)])/2, 
         col = "tomato", lab = expression("SES " * D[pw]), xlab = "Mean SES MPD")

plot.ses(beta.ses$td1[c(T, F)], (ta.ses$mpd.a[seq(1, 95, by = 2)] + ta.ses$mpd.a[seq(2, 96, by = 2)])/2, 
         col = "tomato", lab = expression("SES D'"[pw]), xlab = expression("Mean SES MPD"[a]))


plot.ses(beta.ses$sd0[c(T, F)], (sa.ses$mpd[seq(1, 95, by = 2)] + sa.ses$mpd[seq(2, 96, by = 2)])/2, 
         col = "skyblue", lab = expression("SES " * D[pw]), xlab = "Mean SES MPD")

plot.ses(beta.ses$sd1[c(T, F)], (sa.ses$mpd.a[seq(1, 95, by = 2)] + sa.ses$mpd.a[seq(2, 96, by = 2)])/2, 
         col = "skyblue", lab = expression("SES D'"[pw]), xlab = expression("Mean SES MPD"[a]))


plot.ses(beta.ses$hd0[c(T, F)], (ha.ses$mpd[seq(1, 95, by = 2)] + ha.ses$mpd[seq(2, 96, by = 2)])/2, 
         col = "forestgreen", lab = expression("SES " * D[pw]), xlab = "Mean SES MPD")

plot.ses(beta.ses$hd1[c(T, F)], (ha.ses$mpd.a[seq(1, 95, by = 2)] + ha.ses$mpd.a[seq(2, 96, by = 2)])/2, 
         col = "forestgreen", lab = expression("SES D'"[pw]), xlab = expression("Mean SES MPD"[a]))
par(op)