load("clean.data/traits.rda")

load("clean.data/tree-abund.rda")
load("clean.data/shrub-abund.rda")
load("clean.data/herb-abund.rda")


#____________________________________________________________________________________________________________________

rownames(tr) <- sub(" ", "_", rownames(tr))

tree.tr <- tr[which(rownames(tr) %in% colnames(ta)), ]
tree.tr[, c(9, 11, 18, 19)] <- NULL
tree.tr[13, 16] <- NA
tree.tr[8, 6] <- NA
ta[, which(!colnames(ta) %in% rownames(tree.tr))] <- NULL

shrub.tr <- tr[which(rownames(tr) %in% colnames(sa)), ]
shrub.tr[, c(9, 11, 18)] <- NULL
shrub.tr[32, 17] <- NA
sa[, which(!colnames(sa) %in% rownames(shrub.tr))] <- NULL

herb.tr <- tr[which(rownames(tr) %in% colnames(ha)), ]
herb.tr[, c(9, 20)] <- NULL
ha[, which(!colnames(ha) %in% rownames(herb.tr))] <- NULL


library(cluster)

tqual <- c(7, 8, 10, 11, 13)
hqual <- c(7, 8, 10, 11, 12, 14, 17)

t.dist <- as.matrix(daisy(tree.tr, metric = "gower"))
t.dist.qt <- as.matrix(daisy(tree.tr[, -tqual], metric = "gower"))
t.dist.ql <- as.matrix(daisy(tree.tr[, tqual], metric = "gower"))

s.dist <- as.matrix(daisy(shrub.tr, metric = "gower"))
s.dist.qt <- as.matrix(daisy(shrub.tr[, -tqual], metric = "gower"))
s.dist.ql <- as.matrix(daisy(shrub.tr[, tqual], metric = "gower"))

h.dist <- as.matrix(daisy(herb.tr, metric = "gower"))
h.dist.qt <- as.matrix(daisy(herb.tr[, -hqual], metric = "gower"))
h.dist.ql <- as.matrix(daisy(herb.tr[, hqual], metric = "gower"))


save(t.dist, t.dist.qt, t.dist.ql, file = "clean.data/tree-traits.rda")
save(s.dist, s.dist.qt, s.dist.ql, file = "clean.data/shrub-traits.rda")
save(h.dist, h.dist.qt, h.dist.ql, file = "clean.data/herb-traits.rda")

save(ta, sa, ha, file = "clean.data/func-abund.rda")


#___________________________________________________________________________________________________

library(picante)

t.ses <- as.data.frame(matrix(, nrow = 96, ncol = 0))

ses.mpd <- ses.mpd(ta, t.dist, null.model = "independentswap", runs = 999, iterations = 1000)
t.ses$mpd.z <- ses.mpd$mpd.obs.z
t.ses$mpd.p <- ses.mpd$mpd.obs.p

ses.mpd.a <- ses.mpd(ta, t.dist, null.model = "independentswap", runs = 999, iterations = 1000, abundance.weighted = T)
t.ses$mpd.a.z <- ses.mpd.a$mpd.obs.z
t.ses$mpd.a.p <- ses.mpd.a$mpd.obs.p


t.ses.qt <- as.data.frame(matrix(, nrow = 96, ncol = 0))

ses.mpd <- ses.mpd(ta, t.dist.qt, null.model = "independentswap", runs = 999, iterations = 1000)
t.ses.qt$mpd.z <- ses.mpd$mpd.obs.z
t.ses.qt$mpd.p <- ses.mpd$mpd.obs.p

ses.mpd.a <- ses.mpd(ta, t.dist.qt, null.model = "independentswap", runs = 999, iterations = 1000, abundance.weighted = T)
t.ses.qt$mpd.a.z <- ses.mpd.a$mpd.obs.z
t.ses.qt$mpd.a.p <- ses.mpd.a$mpd.obs.p


t.ses.ql <- as.data.frame(matrix(, nrow = 96, ncol = 0))

ses.mpd <- ses.mpd(ta, t.dist.ql, null.model = "independentswap", runs = 999, iterations = 1000)
t.ses.ql$mpd.z <- ses.mpd$mpd.obs.z
t.ses.ql$mpd.p <- ses.mpd$mpd.obs.p

ses.mpd.a <- ses.mpd(ta, t.dist.ql, null.model = "independentswap", runs = 999, iterations = 1000, abundance.weighted = T)
t.ses.ql$mpd.a.z <- ses.mpd.a$mpd.obs.z
t.ses.ql$mpd.a.p <- ses.mpd.a$mpd.obs.p


save(t.ses, t.ses.qt, t.ses.ql, file = "clean.data/tree-trait-ses.rda")

#__________________________________________________________________________________


s.ses <- as.data.frame(matrix(, nrow = 96, ncol = 0))

ses.mpd <- ses.mpd(sa, s.dist, null.model = "independentswap", runs = 999, iterations = 1000)
s.ses$mpd.z <- ses.mpd$mpd.obs.z
s.ses$mpd.p <- ses.mpd$mpd.obs.p

ses.mpd.a <- ses.mpd(sa, s.dist, null.model = "independentswap", runs = 999, iterations = 1000, abundance.weighted = T)
s.ses$mpd.a.z <- ses.mpd.a$mpd.obs.z
s.ses$mpd.a.p <- ses.mpd.a$mpd.obs.p


s.ses.qt <- as.data.frame(matrix(, nrow = 96, ncol = 0))

ses.mpd <- ses.mpd(sa, s.dist.qt, null.model = "independentswap", runs = 999, iterations = 1000)
s.ses.qt$mpd.z <- ses.mpd$mpd.obs.z
s.ses.qt$mpd.p <- ses.mpd$mpd.obs.p

ses.mpd.a <- ses.mpd(sa, s.dist.qt, null.model = "independentswap", runs = 999, iterations = 1000, abundance.weighted = T)
s.ses.qt$mpd.a.z <- ses.mpd.a$mpd.obs.z
s.ses.qt$mpd.a.p <- ses.mpd.a$mpd.obs.p


s.ses.ql <- as.data.frame(matrix(, nrow = 96, ncol = 0))

ses.mpd <- ses.mpd(sa, s.dist.ql, null.model = "independentswap", runs = 999, iterations = 1000)
s.ses.ql$mpd.z <- ses.mpd$mpd.obs.z
s.ses.ql$mpd.p <- ses.mpd$mpd.obs.p

ses.mpd.a <- ses.mpd(sa, s.dist.ql, null.model = "independentswap", runs = 999, iterations = 1000, abundance.weighted = T)
s.ses.ql$mpd.a.z <- ses.mpd.a$mpd.obs.z
s.ses.ql$mpd.a.p <- ses.mpd.a$mpd.obs.p


save(s.ses, s.ses.qt, s.ses.ql, file = "clean.data/shrub-trait-ses.rda")

#__________________________________________________________________________________


h.ses <- as.data.frame(matrix(, nrow = 96, ncol = 0))

ses.mpd <- ses.mpd(ha, h.dist, null.model = "independentswap", runs = 999, iterations = 1000)
h.ses$mpd.z <- ses.mpd$mpd.obs.z
h.ses$mpd.p <- ses.mpd$mpd.obs.p

ses.mpd.a <- ses.mpd(ha, h.dist, null.model = "independentswap", runs = 999, iterations = 1000, abundance.weighted = T)
h.ses$mpd.a.z <- ses.mpd.a$mpd.obs.z
h.ses$mpd.a.p <- ses.mpd.a$mpd.obs.p


h.ses.qt <- as.data.frame(matrix(, nrow = 96, ncol = 0))

ses.mpd <- ses.mpd(ha, h.dist.qt, null.model = "independentswap", runs = 999, iterations = 1000)
h.ses.qt$mpd.z <- ses.mpd$mpd.obs.z
h.ses.qt$mpd.p <- ses.mpd$mpd.obs.p

ses.mpd.a <- ses.mpd(ha, h.dist.qt, null.model = "independentswap", runs = 999, iterations = 1000, abundance.weighted = T)
h.ses.qt$mpd.a.z <- ses.mpd.a$mpd.obs.z
h.ses.qt$mpd.a.p <- ses.mpd.a$mpd.obs.p


h.ses.ql <- as.data.frame(matrix(, nrow = 96, ncol = 0))

ses.mpd <- ses.mpd(ha, h.dist.ql, null.model = "independentswap", runs = 999, iterations = 1000)
h.ses.ql$mpd.z <- ses.mpd$mpd.obs.z
h.ses.ql$mpd.p <- ses.mpd$mpd.obs.p

ses.mpd.a <- ses.mpd(ha, h.dist.ql, null.model = "independentswap", runs = 999, iterations = 1000, abundance.weighted = T)
h.ses.ql$mpd.a.z <- ses.mpd.a$mpd.obs.z
h.ses.ql$mpd.a.p <- ses.mpd.a$mpd.obs.p


save(h.ses, h.ses.qt, h.ses.ql, file = "clean.data/herb-trait-ses.rda")



