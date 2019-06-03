
library(V.PhyloMaker)
library(dplyr)
library(readxl)


scbi <- read.csv2("raw.data/SCBI_initial_woody_stem_census_2012.csv", sep = ",")

scbi.sp.list <- unique(scbi$Latin)
scbi.sp.list <- sub(" ", "_", scbi.sp.list)
scbi.g.list <- sapply(strsplit(scbi.sp.list, "_"), function(x) paste(x[1]))

id <- which(scbi.sp.list %in% GBOTB.extended$tip.label)
scbi.sp.list[-id]

scbi.sp.list <- scbi.sp.list[-68] #which is "Unidentified unknown"

gf <- tips.info[, c(3,4)] %>% group_by(family, genus) %>% slice(1) %>% ungroup()
scbi.f.list <- c()
for (i in 1:length(scbi.g.list)) {
  scbi.f.list <- c(scbi.f.list, gf$family[which(gf$genus %in% scbi.g.list[i])])
}

scbi.df <- data.frame(species = scbi.sp.list, genus = scbi.g.list, family = scbi.f.list)
scbi.df[, c("species.relative", "genus.relative")] <- NA

rm(id, scbi.sp.list, scbi.g.list, scbi.f.list)

#_____________________________________________________________

load("raw.data/meadows-2018.rda")
rm(meta)

md.sp.list <- colnames(dat)
md.g.list <- sapply(strsplit(md.sp.list, "_"), function(x) paste(x[1]))

id2 <- which(md.sp.list %in% GBOTB.extended$tip.label)
md.sp.list[-id2]
md.sp.list[112] <- "Silene_flos-cuculi"

md.f.list <- c()
for (i in 1:length(md.g.list)) {
  if (md.g.list[i] %in% gf$genus){
    md.f.list <- c(md.f.list, gf$family[which(gf$genus %in% md.g.list[i])])
  } else {
    md.f.list <- c(md.f.list, NA)
  }
}

md.f.list[105:107] <- "Polygonaceae"

md.df <- data.frame(species = md.sp.list, genus = md.g.list, family = md.f.list)
md.df[, c("species.relative", "genus.relative")] <- NA

rm(id2, md.sp.list, md.g.list, md.f.list)

#_____________________________________________________________

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

scbi.tree <- phylo.maker(sp.list = scbi.df, scenarios="S1")
md.tree <- phylo.maker(sp.list = md.df, scenarios="S1")
r.tree <- phylo.maker(sp.list = r.df, scenarios="S1")
rh.tree <- phylo.maker(sp.list = rh.df, scenarios="S1")

scbi.tree <- scbi.tree$scenario.1
md.tree <- md.tree$scenario.1
r.tree <- r.tree$scenario.1
rh.tree <- rh.tree$scenario.1

save(scbi.tree, file = "clean.data/scbi-phylo.rda")
save(md.tree, file = "clean.data/meadows-phylo.rda")
save(r.tree, file = "clean.data/russian-forest-wood-phylo.rda")
save(rh.tree, file = "clean.data/russian-forest-herb-phylo.rda")

save(ta, file = "clean.data/russian-forest-tree-abund.rda")
save(sa, file = "clean.data/russian-forest-shrub-abund.rda")
save(ha, file = "clean.data/russian-forest-herb-abund.rda")
