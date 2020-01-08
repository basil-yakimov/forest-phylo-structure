library(picante)

load("clean.data/scaling.rda")
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
s.dist <- as.matrix(daisy(shrub.tr, metric = "gower"))
h.dist <- as.matrix(daisy(herb.tr, metric = "gower"))



#___________________________________________________________________________________________________



ta_sc_ses <- sa_sc_ses <- ha_sc_ses <- vector(mode = "list", length = 20)

null <- "taxa.labels"

for (jj in 1:20)
{
  print(jj)
  
  ab <- ta_sc[[jj]]
  ab <- ab[, colnames(ta_sc[[jj]]) %in% colnames(t.dist)]
  
  ses <- ses.mpd(ab, t.dist[colnames(t.dist) %in% colnames(ab), colnames(t.dist) %in% colnames(ab)], null.model = null)
  ses.a <- ses.mpd(ab, t.dist[colnames(t.dist) %in% colnames(ab), colnames(t.dist) %in% colnames(ab)], null.model = null,
                   abundance.weighted = T)
  ta_sc_ses[[jj]] <- cbind(z.tl = ses$mpd.obs.z, p.tl = ses$mpd.obs.p,
                           z.a.tl = ses.a$mpd.obs.z, p.a.tl = ses.a$mpd.obs.p)
  
  ab <- sa_sc[[jj]]
  ab <- ab[, colnames(sa_sc[[jj]]) %in% colnames(s.dist)]
  
  ses <- ses.mpd(ab, s.dist[colnames(s.dist) %in% colnames(ab), colnames(s.dist) %in% colnames(ab)], null.model = null)
  ses.a <- ses.mpd(ab, s.dist[colnames(s.dist) %in% colnames(ab), colnames(s.dist) %in% colnames(ab)], null.model = null,
                   abundance.weighted = T)
  sa_sc_ses[[jj]] <- cbind(z.tl = ses$mpd.obs.z, p.tl = ses$mpd.obs.p,
                           z.a.tl = ses.a$mpd.obs.z, p.a.tl = ses.a$mpd.obs.p)
  
  ab <- ha_sc[[jj]]
  ab <- ab[, colnames(ha_sc[[jj]]) %in% colnames(h.dist)]
  
  ses <- ses.mpd(ab, h.dist[colnames(h.dist) %in% colnames(ab), colnames(h.dist) %in% colnames(ab)], null.model = null)
  ses.a <- ses.mpd(ab, h.dist[colnames(h.dist) %in% colnames(ab), colnames(h.dist) %in% colnames(ab)], null.model = null,
                   abundance.weighted = T)
  ha_sc_ses[[jj]] <- cbind(z.tl = ses$mpd.obs.z, p.tl = ses$mpd.obs.p,
                           z.a.tl = ses.a$mpd.obs.z, p.a.tl = ses.a$mpd.obs.p)
}  


save(ta_sc_ses, sa_sc_ses, ha_sc_ses, file = "clean.data/functional-scaling-ses.rda")
