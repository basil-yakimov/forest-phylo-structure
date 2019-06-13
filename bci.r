library(V.PhyloMaker)
library(tidyverse)
library(taRifx)

bci <- read.table("raw.data/bci5.txt", header = T)
tax <- read.csv2("raw.data/Condit_FullBCITaxa.csv")
bcispp <- read.table("raw.data/bcispp.txt", header = T)



bcispp$spp <- paste0(bcispp$genus, " ", bcispp$species)
tax2 <- tax[tax$Status == "Obsolete", ]

for (i in 1:nrow(tax2)){
  if (tax2$Latin[i] %in% bcispp$spp){
    id <- which(bcispp$spp %in% tax2$Latin[i])
    cur <- as.character(tax2$Current[i])
    bcispp$spp[id] <- cur
  }
}

bcispp <- remove.factors(bcispp)
for (i in 1:nrow(tax)){
  if (tax$Current[i] %in% bcispp$spp){
    id <- which(bcispp$spp %in% tax$Current[i])
    f <- as.character(tax$Family[i])
    bcispp$family[id] <- f
  }
}

bcispp$spp <- sapply(strsplit(bcispp$spp, " "), function(x) paste(x[1], x[2]))
bcispp$spp <- sapply(strsplit(bcispp$spp, "_"), function(x) paste(x[1]))
bcispp$spp <- sub(" ", "_", bcispp$spp)
bcispp$family[grep("Fabaceae*", bcispp$family)] <- "Fabaceae"
  
id <- which(bcispp$spp %in% GBOTB.extended$tip.label)
bcispp$spp[-id]
bcispp <- bcispp[-300, ]


for (i in 1:nrow(bcispp)){
  if (bcispp$genus[i] %in% tips.info$genus){
    id <- grep(bcispp$genus[i], tips.info$genus)[1]
    if (!(bcispp$family[i] == tips.info$family[id])){
      bcispp$family[i] <- tips.info$family[id]
    }
  } 
}

bcispp$family[292] <- "Cannabaceae"

sp.list <- data.frame(species = bcispp$spp, genus = bcispp$genus, family = bcispp$family)
bci.tree <- phylo.maker(sp.list = sp.list, scenarios="S1")
#Some taxonomic classifications are not consistent between sp.list and tree

bci.tree <- bci.tree$scenario.1
save(bci.tree, file = "clean.data/bci-phylo.rda")

