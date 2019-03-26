library(vegan)
library(ape)
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

plot(hgt, specnumber(ta), pch = 21, bg = "tomato", 
     xlab = "Высота, м", ylab = "S")

plot(hgt, specnumber(sa), pch = 21, bg = "skyblue", 
     xlab = "Высота, м", ylab = "S")

plot(hgt, specnumber(ha), pch = 21, bg = "forestgreen", 
     xlab = "Высота, м", ylab = "S")

#---#

plot(hgt, exp(diversity(ta)), pch = 21, bg = "tomato", 
     xlab = "Высота, м", ylab = "S")

plot(hgt, exp(diversity(sa)), pch = 21, bg = "skyblue", 
     xlab = "Высота, м", ylab = "S")

plot(hgt, exp(diversity(ha)), pch = 21, bg = "forestgreen", 
     xlab = "Высота, м", ylab = "S")

#---#

plot(hgt, pd(ta, wood.tree)$PD, pch = 21, bg = "tomato", 
     xlab = "Высота, м", ylab = "PD")

plot(hgt, pd(sa, wood.tree)$PD, pch = 21, bg = "skyblue", 
     xlab = "Высота, м", ylab = "S")

plot(hgt, pd(ha, herb.tree)$PD, pch = 21, bg = "forestgreen", 
     xlab = "Высота, м", ylab = "S")



