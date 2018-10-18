library(ape)
phylo <- read.tree("Vascular_Plants_rooted.dated.tre")

library(readxl)
sp1 <- read_excel("Forest_data.xlsx", sheet = 4, range = "A2:A21", col_names = F)[[1]]
sp2 <- read_excel("Forest_data.xlsx", sheet = 5, range = "A2:A42", col_names = F)[[1]]

tree.list <- sapply(strsplit(sp1, " "), function(x) paste(x[1], x[2]))
tree.list <- sub(" ", "_", tree.list)

shrub.list <- sapply(strsplit(sp2, " "), function(x) paste(x[1], x[2]))
shrub.list <- sub(" ", "_", shrub.list)


ta <- read_excel("Forest_data.xlsx", sheet = 4, range = "B2:CS21", col_names = F)
ta <- t(ta)
colnames(ta) <- tree.list
ta <- ta[, colSums(ta) > 0]
ta <- data.frame(ta)

sa <- read_excel("Forest_data.xlsx", sheet = 5, range = "B2:CS42", col_names = F)
sa <- t(sa)
colnames(sa) <- shrub.list
sa <- sa[, colSums(sa) > 0]
sa <- data.frame(sa)

#---#

sp.list <- sort(unique(c(colnames(ta), colnames(sa))))

genera <- sapply(strsplit(sp.list, "_"), function(x) paste(x[1]))
table(genera)

#---#

source("addition-tools.r")

id <- which(sp.list %in% phylo$tip.label)

sp.list[id]
sp.list[-id]

grep("Zabelia*", phylo$tip.label, value = T)
phylo$tip.label[phylo$tip.label == "Zabelia_biflora"] <- "Abelia_biflora"   # Zanne mistake

grep("Betula*", phylo$tip.label, value = T)
phylo$tip.label[phylo$tip.label == "Betula_davurica"] <- "Betula_dahurica"   # Zanne mistake

grep("Carpinus*", phylo$tip.label, value = T)
phylo$tip.label[phylo$tip.label == "Carpinus_turczaninovii"] <- "Carpinus_turczaninowii"   # Zanne mistake

grep("Deutzia*", phylo$tip.label, value = T)
phylo$tip.label[phylo$tip.label == "Deutzia_scabra"] <- "Deutzia_grandiflora"   # two closest Deutsia species
phylo$tip.label[phylo$tip.label == "Deutzia_rubens"] <- "Deutzia_parviflora"    # two closest Deutsia species

t(t(sort(grep("Populus*", phylo$tip.label, value = T))))
phylo$tip.label[phylo$tip.label == "Populus_trichocarpa"] <- "Populus_cathayana"   # one of species in the same subsection Tacamahaca; see the phylogeny at doi: 10.3732/ajb.91.9.1398 
sp.list[sp.list == "Populus_davidiana"] <- "Populus_tremula"   # Yuxin synonim: subspecies of Populus tremula
colnames(sa)[colnames(sa) == "Populus_davidiana"] <- "Populus_tremula"
colnames(ta)[colnames(ta) == "Populus_davidiana"] <- "Populus_tremula"

t(t(sort(grep("Quercus*", phylo$tip.label, value = T))))
phylo$tip.label[phylo$tip.label == "Quercus_mongolica"] <- "Quercus_liaotungensis"   # one of the closest species

t(t(sort(grep("Rhamnus*", phylo$tip.label, value = T))))
phylo$tip.label[phylo$tip.label == "Rhamnus_saxatilis"] <- "Rhamnus_utilis"   # a species close to R.davurica; see also DOI: 10.2307/4135616 

t(t(sort(grep("Rubus*", phylo$tip.label, value = T))))
phylo$tip.label[phylo$tip.label == "Rubus_idaeus"] <- "Rubus_Sp." # any Rubus

t(t(sort(grep("Spiraea*", phylo$tip.label, value = T))))
phylo <- my.add(phylo, "Spiraea_dasyantha", tol = 1e-6) # add polytomy to genus

id1 <- grep("Spiraea*", phylo$tip.label)
tree <- drop.tip(phylo, phylo$tip.label[-id1])
plot(tree, cex = 0.7)

t(t(sort(grep("Syringa*", phylo$tip.label, value = T))))
sp.list <- sp.list[sp.list != "Syringa_pekinensis"]      # Yuxin synonim: subspecies of Syringa reticulata
sa$Syringa_reticulata <- sa$Syringa_reticulata + sa$Syringa_pekinensis
sa$Syringa_pekinensis <- NULL

t(t(sort(grep("Tilia*", phylo$tip.label, value = T))))
phylo$tip.label[phylo$tip.label == "Tilia_cordata"] <- "Tilia_amurensis"         # species close to T.cordata, see DOI: 10.1111/j.1095-8339.2008.00891.x

t(t(sort(grep("Ulmus*", phylo$tip.label, value = T))))
phylo$tip.label[phylo$tip.label == "Ulmus_glabra"] <- "Ulmus_japonica"  # species from the same section, see wiki and DOI: 10.2307/2419779 
phylo <- my.add(phylo, "Ulmus_laciniata", tol = 1e-6) # add polytomy to genus

id1 <- grep("Ulmus*", phylo$tip.label)
tree <- drop.tip(phylo, phylo$tip.label[-id1])
plot(tree, cex = 0.7)

t(t(sort(grep("Viburnum*", phylo$tip.label, value = T))))
t(t(sort(grep("Lonicera*", phylo$tip.label, value = T))))
phylo$tip.label[phylo$tip.label == "Viburnum_opulus"] <- "Viburnum_mongolicum" # any Viburnum

#-------#

id <- which(sp.list %in% phylo$tip.label)

sp.list[id]
sp.list[-id]

id <- which(phylo$tip.label %in% sp.list)
tree <- drop.tip(phylo, phylo$tip.label[-id])
plot(tree, cex = 0.7)
axisPhylo()

is.ultrametric(tree, tol = 1e-8)

#-------#

rm(phylo, genera, shrub.list, tree.list, sp1, sp2, sp.list, id, id1)
save.image("tree-shrub-phylo.RData")