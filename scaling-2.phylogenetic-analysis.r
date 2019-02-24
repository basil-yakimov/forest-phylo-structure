library(ape)
library(picante)

load("clean.data/wood-phylo.rda")
load("clean.data/herb-phylo.rda")

wood.mat <- cophenetic(wood.tree)
herb.mat <- cophenetic(herb.tree)

load("clean.data/scaling-abund.rda")


#________________________________________________

c <- vector(mode = "list", length = 3)
for (ii in 1:3){
  
  g <- vector(mode = "list", length = 7)
  
  for (jj in 1:7){
    df <- out1_tot[[ii]][[jj]]
    
    if (ii == 3){
      ses <- ses.mpd(df, herb.mat, null.model = "richness", runs = 999, iterations = 1000)
    } else {
      ses <- ses.mpd(df, wood.mat, null.model = "richness", runs = 999, iterations = 1000)
    }
    m <- cbind(ses$mpd.obs.z, ses$mpd.obs.p)
    colnames(m) <- c("mpd.obs.z", "mpd.obs.p")
    g[[jj]] <- m
    
  }
  
  c[[ii]] <- g
}

#________________________________________________

c2 <- vector(mode = "list", length = 3)
for (ii in 1:3){
  
  g <- vector(mode = "list", length = 95)
  
  for (jj in 1:7){
    df <- out1_tot[[ii]][[jj]]
    
    if (ii == 3){
      ses <- ses.mpd(df, herb.mat, null.model = "richness", runs = 999, iterations = 1000)
    } else {
      ses <- ses.mpd(df, wood.mat, null.model = "richness", runs = 999, iterations = 1000)
    }
    m <- cbind(ses$mpd.obs.z, ses$mpd.obs.p)
    colnames(m) <- c("mpd.obs.z", "mpd.obs.p")
    g[[jj]] <- m
    
  }
  
  c2[[ii]] <- g
}

save(c, c2, file = "clean.data/scaling-phylo-ses.rda")
