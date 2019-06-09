library(V.PhyloMaker)
library(dplyr)
library(readxl)



#wood layers

sp1 <- read_excel("raw.data/Forest_data.xlsx", sheet = 1, range = "A2:A12", col_names = F)[[1]]
sp2 <- read_excel("raw.data/Forest_data.xlsx", sheet = 2, range = "A2:A11", col_names = F)[[1]]

tree.list <- sapply(strsplit(sp1, " "), function(x) paste(x[1], x[2]))
tree.list <- sub(" ", "_", tree.list)
tree.list[5] <- "Prunus_cerasifera" #synonym

shrub.list <- sapply(strsplit(sp2, " "), function(x) paste(x[1], x[2]))
shrub.list <- sub(" ", "_", shrub.list)
shrub.list[6] <- "Prunus_cerasifera"

rm(sp1, sp2)


ta <- read_excel("raw.data/Forest_data.xlsx", sheet = 1, range = "B2:CW12", col_names = F)
ta <- t(ta)
colnames(ta) <- tree.list
ta <- ta[, colSums(ta) > 0]
ta <- data.frame(ta)

sa <- read_excel("raw.data/Forest_data.xlsx", sheet = 2, range = "B2:CW11", col_names = F)
sa <- t(sa)
colnames(sa) <- shrub.list
sa <- sa[, colSums(sa) > 0]
sa <- data.frame(sa)

rm(tree.list, shrub.list)


r.sp.list <- sort(unique(c(colnames(ta), colnames(sa))))
r.g.list <- sapply(strsplit(r.sp.list, "_"), function(x) paste(x[1]))


id3 <- which(r.sp.list %in% GBOTB.extended$tip.label)
r.sp.list[-id3]

gf <- tips.info[, c(3,4)] %>% group_by(family, genus) %>% slice(1) %>% ungroup()
r.f.list <- c()
for (i in 1:length(r.g.list)) {
  if (r.g.list[i] %in% gf$genus){
    r.f.list <- c(r.f.list, gf$family[which(gf$genus %in% r.g.list[i])])
  } else {
    r.f.list <- c(r.f.list, NA)
  }
}

r.df <- data.frame(species = r.sp.list, genus = r.g.list, family = r.f.list)
r.df[, c("species.relative", "genus.relative")] <- NA

rm(id3, r.sp.list, r.g.list, r.f.list)

#_____________________________________________________________

#herb layer

sp3 <- read_excel("raw.data/Forest_data.xlsx", sheet = 3, range = "A2:A43", col_names = F)[[1]]

herb.list <- sapply(strsplit(sp3, " "), function(x) paste(x[1], x[2]))
herb.list <- sub(" ", "_", herb.list)

ha <- read_excel("raw.data/Forest_data.xlsx", sheet = 3, range = "B2:CW43", col_names = F)
ha <- t(ha)
colnames(ha) <- herb.list
ha <- ha[, colSums(ha) > 0]
ha <- data.frame(ha)

herb.g.list <- sapply(strsplit(herb.list, "_"), function(x) paste(x[1]))

id4 <- which(herb.list %in% GBOTB.extended$tip.label)
herb.list[-id4]

herb.f.list <- c()
for (i in 1:length(herb.g.list)) {
  if (herb.g.list[i] %in% gf$genus){
    herb.f.list <- c(herb.f.list, gf$family[which(gf$genus %in% herb.g.list[i])])
  } else {
    herb.f.list <- c(herb.f.list, NA)
  }
}

rh.df <- data.frame(species = herb.list, genus = herb.g.list, family = herb.f.list)
rh.df[, c("species.relative", "genus.relative")] <- NA

rm(id4, sp3, herb.list, herb.g.list, herb.f.list)

#_____________________________________________________________

r.tree <- phylo.maker(sp.list = r.df, scenarios="S1")
rh.tree <- phylo.maker(sp.list = rh.df, scenarios="S1")

r.tree <- r.tree$scenario.1
rh.tree <- rh.tree$scenario.1

rh.tree$tip.label[grep("Athyrium_filix-femina", rh.tree$tip.label)] <- "Athyrium_filix.femina"
rh.tree$tip.label[grep("Dryopteris_filix-mas", rh.tree$tip.label)] <- "Dryopteris_filix.mas"

save(r.tree, file = "clean.data/russian-forest-wood-phylo.rda")
save(rh.tree, file = "clean.data/russian-forest-herb-phylo.rda")

#_____________________________________________________________

hgt <- read.table("raw.data/hgt.txt")[[1]]

sc <- 1:20

hgt_sc <- ta_sc <- sa_sc <- ha_sc <- vector(mode = "list", length = length(sc)) 

for (jj in 1:length(sc))
{
  num <- 96 %/% sc[jj]
  
  dt <- matrix(0, nrow = num, ncol = ncol(ta))
  ds <- matrix(0, nrow = num, ncol = ncol(sa))
  dh <- matrix(0, nrow = num, ncol = ncol(ha))
  dhgt <- rep(NA, num)
  
  colnames(dt) <- colnames(ta)
  colnames(ds) <- colnames(sa)
  colnames(dh) <- colnames(ha)
  
  for (ii in 1:num)
  {
    dt[ii, ] <- colSums(ta[((ii-1)*sc[jj]+1):(ii*sc[jj]), ])
    ds[ii, ] <- colSums(sa[((ii-1)*sc[jj]+1):(ii*sc[jj]), ])
    dh[ii, ] <- colSums(ha[((ii-1)*sc[jj]+1):(ii*sc[jj]), ])
    dhgt[ii] <- mean(hgt[((ii-1)*sc[jj]+1):(ii*sc[jj])])
  }
  
  ta_sc[[jj]] <- dt
  sa_sc[[jj]] <- ds
  ha_sc[[jj]] <- dh
  hgt_sc[[jj]] <- dhgt
}

save(hgt_sc, ta_sc, sa_sc, ha_sc, file = "clean.data/russian-forest-scaling.rda")

#_____________________________________________________________

library(ape)
library(picante)

load("clean.data/russian-forest-wood-phylo.rda")
load("clean.data/russian-forest-herb-phylo.rda")

load("clean.data/russian-forest-scaling.rda")

wood.mat <- cophenetic(r.tree)
herb.mat <- cophenetic(rh.tree)

ta_sc_ses <- sa_sc_ses <- ha_sc_ses <- vector(mode = "list", length = 20)

null <- "taxa.labels"

for (jj in 1:20)
{
  print(jj)
  
  ses <- ses.mpd(ta_sc[[jj]], wood.mat, null.model = null)
  ses.a <- ses.mpd(ta_sc[[jj]], wood.mat, null.model = null, abundance.weighted = T)
  ta_sc_ses[[jj]] <- cbind(z.tl = ses$mpd.obs.z, p.tl = ses$mpd.obs.p,
                           z.a.tl = ses.a$mpd.obs.z, p.a.tl = ses.a$mpd.obs.p)
  
  ses <- ses.mpd(sa_sc[[jj]], wood.mat, null.model = null)
  ses.a <- ses.mpd(sa_sc[[jj]], wood.mat, null.model = null, abundance.weighted = T)
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
  
  ses <- ses.mpd(sa_sc[[jj]], wood.mat, null.model = null)
  ses.a <- ses.mpd(sa_sc[[jj]], wood.mat, null.model = null, abundance.weighted = T)
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
  
  ses <- ses.mpd(sa_sc[[jj]], wood.mat, null.model = null)
  ses.a <- ses.mpd(sa_sc[[jj]], wood.mat, null.model = null, abundance.weighted = T)
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
  
  ses <- ses.mpd(sa_sc[[jj]], wood.mat, null.model = null)
  ses.a <- ses.mpd(sa_sc[[jj]], wood.mat, null.model = null, abundance.weighted = T)
  sa_sc_ses[[jj]] <- cbind(sa_sc_ses[[jj]], cbind(z.ts = ses$mpd.obs.z, p.ts = ses$mpd.obs.p,
                                                  z.a.ts = ses.a$mpd.obs.z, p.a.ts = ses.a$mpd.obs.p))
  
  ses <- ses.mpd(ha_sc[[jj]], herb.mat, null.model = null)
  ses.a <- ses.mpd(ha_sc[[jj]], herb.mat, null.model = null, abundance.weighted = T)
  ha_sc_ses[[jj]] <- cbind(ha_sc_ses[[jj]], cbind(z.ts = ses$mpd.obs.z, p.ts = ses$mpd.obs.p,
                                                  z.a.ts = ses.a$mpd.obs.z, p.a.ts = ses.a$mpd.obs.p))
}  

save(ta_sc_ses, sa_sc_ses, ha_sc_ses, file = "clean.data/russian-forest-scaling-ses.rda")

#_____________________________________________________________

load("clean.data/russian-forest-scaling.rda")
load("clean.data/russian-forest-scaling-ses.rda")
sc <- 1:20


png("figures/rf-mean-nri-sc.png", 4000, 4000, pointsize = 75)

op <- par(mfrow = c(4,3), mar = c(3, 3, 2, 0.5))

meanNRI <- sapply(ta_sc_ses, function(x) -mean(x[, "z.tl"], na.rm = T))
plot(sc, meanNRI, pch = 21, bg = "tomato", main = "Tree: taxa labels")
model <- lm(meanNRI ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model)

meanNRI <- sapply(sa_sc_ses, function(x) -mean(x[, "z.tl"], na.rm = T))
plot(sc, meanNRI, pch = 21, bg = "skyblue", main = "Shrub: taxa labels")
model <- lm(meanNRI ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model)

meanNRI <- sapply(ha_sc_ses, function(x) -mean(x[, "z.tl"], na.rm = T))
plot(sc, meanNRI, pch = 21, bg = "forestgreen", main = "Herb: taxa labels")
model <- lm(meanNRI ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model)

#---#

meanNRI <- sapply(ta_sc_ses, function(x) -mean(x[, "z.r"], na.rm = T))
plot(sc, meanNRI, pch = 21, bg = "tomato", main = "Tree: richness")
model <- lm(meanNRI ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model)

meanNRI <- sapply(sa_sc_ses, function(x) -mean(x[, "z.r"], na.rm = T))
plot(sc, meanNRI, pch = 21, bg = "skyblue", main = "Shrub: richness")
model <- lm(meanNRI ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model)

meanNRI <- sapply(ha_sc_ses, function(x) -mean(x[, "z.r"], na.rm = T))
plot(sc, meanNRI, pch = 21, bg = "forestgreen", main = "Herb: richness")
model <- lm(meanNRI ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model)

#---#

meanNRI <- sapply(ta_sc_ses, function(x) -mean(x[, "z.is"], na.ism = T))
plot(sc, meanNRI, pch = 21, bg = "tomato", main = "Tree: ind. swap")
model <- lm(meanNRI ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model)

meanNRI <- sapply(sa_sc_ses, function(x) -mean(x[, "z.is"], na.ism = T))
plot(sc, meanNRI, pch = 21, bg = "skyblue", main = "Shrub: ind. swap")
model <- lm(meanNRI ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model)

meanNRI <- sapply(ha_sc_ses, function(x) -mean(x[, "z.is"], na.ism = T))
plot(sc, meanNRI, pch = 21, bg = "forestgreen", main = "Herb: ind. swap")
model <- lm(meanNRI ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model)

#---#

meanNRI <- sapply(ta_sc_ses, function(x) -mean(x[, "z.ts"], na.ism = T))
plot(sc, meanNRI, pch = 21, bg = "tomato", main = "Tree: trial swap")
model <- lm(meanNRI ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model)

meanNRI <- sapply(sa_sc_ses, function(x) -mean(x[, "z.ts"], na.ism = T))

plot(sc, meanNRI, pch = 21, bg = "skyblue", main = "Shrub: trial swap")
model <- lm(meanNRI ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model)

meanNRI <- sapply(ha_sc_ses, function(x) -mean(x[, "z.ts"], na.ism = T))

plot(sc, meanNRI, pch = 21, bg = "forestgreen", main = "Herb: trial swap")
model <- lm(meanNRI ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model)

dev.off()



#-----------------------------------------------------------------------------#

png("figures/rf-mean-nri-a-sc.png", 4000, 4000, pointsize = 75)

op <- par(mfrow = c(4,3), mar = c(3, 3, 2, 0.5))

meanNRI <- sapply(ta_sc_ses, function(x) -mean(x[, "z.a.tl"], na.rm = T))
plot(sc, meanNRI, pch = 21, bg = "tomato", main = "Tree: taxa labels")
model <- lm(meanNRI ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model)

meanNRI <- sapply(sa_sc_ses, function(x) -mean(x[, "z.a.tl"], na.rm = T))
plot(sc, meanNRI, pch = 21, bg = "skyblue", main = "Shrub: taxa labels")
model <- lm(meanNRI ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model)

meanNRI <- sapply(ha_sc_ses, function(x) -mean(x[, "z.a.tl"], na.rm = T))
plot(sc, meanNRI, pch = 21, bg = "forestgreen", main = "Herb: taxa labels")
model <- lm(meanNRI ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model)

#---#

meanNRI <- sapply(ta_sc_ses, function(x) -mean(x[, "z.a.r"], na.rm = T))
plot(sc, meanNRI, pch = 21, bg = "tomato", main = "Tree: richness")
model <- lm(meanNRI ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model)

meanNRI <- sapply(sa_sc_ses, function(x) -mean(x[, "z.a.r"], na.rm = T))
plot(sc, meanNRI, pch = 21, bg = "skyblue", main = "Shrub: richness")
model <- lm(meanNRI ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model)

meanNRI <- sapply(ha_sc_ses, function(x) -mean(x[, "z.a.r"], na.rm = T))
plot(sc, meanNRI, pch = 21, bg = "forestgreen", main = "Herb: richness")
model <- lm(meanNRI ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model)

#---#

meanNRI <- sapply(ta_sc_ses, function(x) -mean(x[, "z.a.is"], na.ism = T))
plot(sc, meanNRI, pch = 21, bg = "tomato", main = "Tree: ind. swap")
model <- lm(meanNRI ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model)

meanNRI <- sapply(sa_sc_ses, function(x) -mean(x[, "z.a.is"], na.ism = T))
plot(sc, meanNRI, pch = 21, bg = "skyblue", main = "Shrub: ind. swap")
model <- lm(meanNRI ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model)

meanNRI <- sapply(ha_sc_ses, function(x) -mean(x[, "z.a.is"], na.ism = T))
plot(sc, meanNRI, pch = 21, bg = "forestgreen", main = "Herb: ind. swap")
model <- lm(meanNRI ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model)

#---#

meanNRI <- sapply(ta_sc_ses, function(x) -mean(x[, "z.a.ts"], na.ism = T))
plot(sc, meanNRI, pch = 21, bg = "tomato", main = "Tree: trial swap")
model <- lm(meanNRI ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model)

meanNRI <- sapply(sa_sc_ses, function(x) -mean(x[, "z.a.ts"], na.ism = T))

plot(sc, meanNRI, pch = 21, bg = "skyblue", main = "Shrub: trial swap")
model <- lm(meanNRI ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model)

meanNRI <- sapply(ha_sc_ses, function(x) -mean(x[, "z.a.ts"], na.ism = T))

plot(sc, meanNRI, pch = 21, bg = "forestgreen", main = "Herb: trial swap")
model <- lm(meanNRI ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model)

dev.off()


