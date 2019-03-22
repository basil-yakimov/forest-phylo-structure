library(ape)
library(picante)

load("clean.data/wood-phylo.rda")
load("clean.data/herb-phylo.rda")

wood.mat <- cophenetic(wood.tree)
herb.mat <- cophenetic(herb.tree)

load("clean.data/scaling.rda")

#________________________________________________

ta_sc_ses <- sa_sc_ses <- ha_sc_ses <- vector(mode = "list", length = 20)

for (jj in 1:20)
{
  ses <- ses.mpd(ta_sc[[jj]], wood.mat, null.model = "richness")
  ses.a <- ses.mpd(ta_sc[[jj]], wood.mat, null.model = "richness", abundance.weighted = T)
  ta_sc_ses[[jj]] <- cbind(z = ses$mpd.obs.z, p = ses$mpd.obs.p,
                           z.a = ses.a$mpd.obs.z, p.a = ses.a$mpd.obs.p)
  
  ses <- ses.mpd(sa_sc[[jj]], wood.mat, null.model = "richness")
  ses.a <- ses.mpd(sa_sc[[jj]], wood.mat, null.model = "richness", abundance.weighted = T)
  sa_sc_ses[[jj]] <- cbind(z = ses$mpd.obs.z, p = ses$mpd.obs.p,
                           z.a = ses.a$mpd.obs.z, p.a = ses.a$mpd.obs.p)
  
  ses <- ses.mpd(ha_sc[[jj]], herb.mat, null.model = "richness")
  ses.a <- ses.mpd(ha_sc[[jj]], herb.mat, null.model = "richness", abundance.weighted = T)
  ha_sc_ses[[jj]] <- cbind(z = ses$mpd.obs.z, p = ses$mpd.obs.p,
                           z.a = ses.a$mpd.obs.z, p.a = ses.a$mpd.obs.p)
}  
  
save(ta_sc_ses, sa_sc_ses, ha_sc_ses, file = "clean.data/scaling-ses.rda")