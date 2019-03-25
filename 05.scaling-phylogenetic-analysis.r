library(ape)
library(picante)

load("clean.data/wood-phylo.rda")
load("clean.data/herb-phylo.rda")

load("clean.data/tree-abund.rda")
load("clean.data/shrub-abund.rda")

wood.mat <- cophenetic(drop.tip(wood.tree, wood.tree$tip.label[!(wood.tree$tip.label %in% colnames(ta))]))
shrub.mat <- cophenetic(drop.tip(wood.tree, wood.tree$tip.label[!(wood.tree$tip.label %in% colnames(sa))]))
herb.mat <- cophenetic(herb.tree)

load("clean.data/scaling.rda")

#________________________________________________

ta_sc_ses <- sa_sc_ses <- ha_sc_ses <- vector(mode = "list", length = 20)

null <- "taxa.labels"

for (jj in 1:20)
{
  print(jj)
  
  ses <- ses.mpd(ta_sc[[jj]], wood.mat, null.model = null)
  ses.a <- ses.mpd(ta_sc[[jj]], wood.mat, null.model = null, abundance.weighted = T)
  ta_sc_ses[[jj]] <- cbind(z.tl = ses$mpd.obs.z, p.tl = ses$mpd.obs.p,
                           z.a.tl = ses.a$mpd.obs.z, p.a.tl = ses.a$mpd.obs.p)
  
  ses <- ses.mpd(sa_sc[[jj]], shrub.mat, null.model = null)
  ses.a <- ses.mpd(sa_sc[[jj]], shrub.mat, null.model = null, abundance.weighted = T)
  sa_sc_ses[[jj]] <- cbind(z.tl = ses$mpd.obs.z, p.tl = ses$mpd.obs.p,
                           z.a.tl = ses.a$mpd.obs.z, p.a.tl = ses.a$mpd.obs.p)
  
  ses <- ses.mpd(ha_sc[[jj]], herb.mat, null.model = null)
  ses.a <- ses.mpd(ha_sc[[jj]], herb.mat, null.model = null, abundance.weighted = T)
  ha_sc_ses[[jj]] <- cbind(z.tl = ses$mpd.obs.z, p.tl = ses$mpd.obs.p,
                           z.a.tl = ses.a$mpd.obs.z, p.a.tl = ses.a$mpd.obs.p)
}  


null <- "richness"

for (jj in 1:20)
{
  print(jj)
  
  ses <- ses.mpd(ta_sc[[jj]], wood.mat, null.model = null)
  ses.a <- ses.mpd(ta_sc[[jj]], wood.mat, null.model = null, abundance.weighted = T)
  ta_sc_ses[[jj]] <- cbind(ta_sc_ses[[jj]], cbind(z.r = ses$mpd.obs.z, p.r = ses$mpd.obs.p,
                           z.a.r = ses.a$mpd.obs.z, p.a.r = ses.a$mpd.obs.p))
  
  ses <- ses.mpd(sa_sc[[jj]], shrub.mat, null.model = null)
  ses.a <- ses.mpd(sa_sc[[jj]], shrub.mat, null.model = null, abundance.weighted = T)
  sa_sc_ses[[jj]] <- cbind(sa_sc_ses[[jj]], cbind(z.r = ses$mpd.obs.z, p.r = ses$mpd.obs.p,
                                                  z.a.r = ses.a$mpd.obs.z, p.a.r = ses.a$mpd.obs.p))
  
  ses <- ses.mpd(ha_sc[[jj]], herb.mat, null.model = null)
  ses.a <- ses.mpd(ha_sc[[jj]], herb.mat, null.model = null, abundance.weighted = T)
  ha_sc_ses[[jj]] <- cbind(ha_sc_ses[[jj]], cbind(z.r = ses$mpd.obs.z, p.r = ses$mpd.obs.p,
                                                  z.a.r = ses.a$mpd.obs.z, p.a.r = ses.a$mpd.obs.p))
}  
  
null <- "independentswap"

for (jj in 1:20)
{
  print(jj)
  
  ses <- ses.mpd(ta_sc[[jj]], wood.mat, null.model = null)
  ses.a <- ses.mpd(ta_sc[[jj]], wood.mat, null.model = null, abundance.weighted = T)
  ta_sc_ses[[jj]] <- cbind(ta_sc_ses[[jj]], cbind(z.is = ses$mpd.obs.z, p.is = ses$mpd.obs.p,
                                                  z.a.is = ses.a$mpd.obs.z, p.a.is = ses.a$mpd.obs.p))
  
  ses <- ses.mpd(sa_sc[[jj]], shrub.mat, null.model = null)
  ses.a <- ses.mpd(sa_sc[[jj]], shrub.mat, null.model = null, abundance.weighted = T)
  sa_sc_ses[[jj]] <- cbind(sa_sc_ses[[jj]], cbind(z.is = ses$mpd.obs.z, p.is = ses$mpd.obs.p,
                                                  z.a.is = ses.a$mpd.obs.z, p.a.is = ses.a$mpd.obs.p))
  
  ses <- ses.mpd(ha_sc[[jj]], herb.mat, null.model = null)
  ses.a <- ses.mpd(ha_sc[[jj]], herb.mat, null.model = null, abundance.weighted = T)
  ha_sc_ses[[jj]] <- cbind(ha_sc_ses[[jj]], cbind(z.is = ses$mpd.obs.z, p.is = ses$mpd.obs.p,
                                                  z.a.is = ses.a$mpd.obs.z, p.a.is = ses.a$mpd.obs.p))
}  

null <- "trialswap"

for (jj in 1:20)
{
  print(jj)
  
  ses <- ses.mpd(ta_sc[[jj]], wood.mat, null.model = null)
  ses.a <- ses.mpd(ta_sc[[jj]], wood.mat, null.model = null, abundance.weighted = T)
  ta_sc_ses[[jj]] <- cbind(ta_sc_ses[[jj]], cbind(z.ts = ses$mpd.obs.z, p.ts = ses$mpd.obs.p,
                                                  z.a.ts = ses.a$mpd.obs.z, p.a.ts = ses.a$mpd.obs.p))
  
  ses <- ses.mpd(sa_sc[[jj]], shrub.mat, null.model = null)
  ses.a <- ses.mpd(sa_sc[[jj]], shrub.mat, null.model = null, abundance.weighted = T)
  sa_sc_ses[[jj]] <- cbind(sa_sc_ses[[jj]], cbind(z.ts = ses$mpd.obs.z, p.ts = ses$mpd.obs.p,
                                                  z.a.ts = ses.a$mpd.obs.z, p.a.ts = ses.a$mpd.obs.p))
  
  ses <- ses.mpd(ha_sc[[jj]], herb.mat, null.model = null)
  ses.a <- ses.mpd(ha_sc[[jj]], herb.mat, null.model = null, abundance.weighted = T)
  ha_sc_ses[[jj]] <- cbind(ha_sc_ses[[jj]], cbind(z.ts = ses$mpd.obs.z, p.ts = ses$mpd.obs.p,
                                                  z.a.ts = ses.a$mpd.obs.z, p.a.ts = ses.a$mpd.obs.p))
}  

save(ta_sc_ses, sa_sc_ses, ha_sc_ses, file = "clean.data/scaling-ses.rda")
