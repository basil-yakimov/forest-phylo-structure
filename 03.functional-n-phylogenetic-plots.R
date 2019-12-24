load("clean.data/tree-trait-ses.rda")
load("clean.data/shrub-trait-ses.rda")
load("clean.data/herb-trait-ses.rda")


op <- par(mfcol = c(3, 2), mar = c(4, 4.1, 1, 1), cex = 0.9)

boxplot(t.ses$mpd.z, s.ses$mpd.z, h.ses$mpd.z, col = c("tomato", "skyblue", "forestgreen"), ylab = "SES MPD", names = c("tree", "shrub", "herb"))
title("all traits")
boxplot(t.ses.ql$mpd.z, s.ses.ql$mpd.z, h.ses.ql$mpd.z, col = c("tomato", "skyblue", "forestgreen"), ylab = "SES MPD", names = c("tree", "shrub", "herb"))
title("qualitative traits")
boxplot(t.ses.qt$mpd.z, s.ses.qt$mpd.z, h.ses.qt$mpd.z, col = c("tomato", "skyblue", "forestgreen"), ylab = "SES MPD", names = c("tree", "shrub", "herb"))
title("quantitative traits")


boxplot(t.ses$mpd.a.z, s.ses$mpd.a.z, h.ses$mpd.a.z, col = c("tomato", "skyblue", "forestgreen"), ylab = expression("SES "*MPD[a]),
        names = c("tree", "shrub", "herb"))
title("all traits")

boxplot(t.ses.ql$mpd.a.z, s.ses.ql$mpd.a.z, h.ses.ql$mpd.a.z, col = c("tomato", "skyblue", "forestgreen"), ylab = expression("SES "*MPD[a]),
        names = c("tree", "shrub", "herb"))
title("qualitative traits")

boxplot(t.ses.qt$mpd.a.z, s.ses.qt$mpd.a.z, h.ses.qt$mpd.a.z, col = c("tomato", "skyblue", "forestgreen"), ylab = expression("SES "*MPD[a]),
        names = c("tree", "shrub", "herb"))
title("quantitative traits")

par(op)

#_____________________________________________________________________________________________________________


load("clean.data/phylo-ses.rda")


op <- par(mfcol = c(1, 2))

boxplot(t.ses.z$independentswap, s.ses.z$independentswap, h.ses.z$independentswap, col = c("tomato", "skyblue", "forestgreen"), ylab = "SES MPD", names = c("tree", "shrub", "herb"))
boxplot(t.ses.z$independentswap.a, s.ses.z$independentswap.a, h.ses.z$independentswap.a, col = c("tomato", "skyblue", "forestgreen"), ylab = expression("SES "*MPD[a]),
        names = c("tree", "shrub", "herb"))

par(op)


