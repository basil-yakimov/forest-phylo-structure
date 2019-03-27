library(vegan)
library(picante)

load("clean.data/wood-phylo.rda")
load("clean.data/herb-phylo.rda")

wood.mat <- cophenetic(wood.tree)
herb.mat <- cophenetic(herb.tree)

load("clean.data/tree-abund.rda")
load("clean.data/shrub-abund.rda")
load("clean.data/herb-abund.rda")

hgt <- read.table("raw.data/hgt.txt")[[1]]

#---#

library(ecomf)

tpd <- PqD(ta, wood.tree, q = 0)

res <- oecosimu(ta, PqD, method = "abuswap_c", nsimul = 9, tree = wood.tree, q = 0)
