
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

scbi.tree <- phylo.maker(sp.list = scbi.df, scenarios="S1")
md.tree <- phylo.maker(sp.list = md.df, scenarios="S1")

scbi.tree <- scbi.tree$scenario.1
md.tree <- md.tree$scenario.1

save(scbi.tree, file = "clean.data/scbi-phylo.rda")
save(md.tree, file = "clean.data/meadows-phylo.rda")
