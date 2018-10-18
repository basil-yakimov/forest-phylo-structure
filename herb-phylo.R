library(ape)
phylo <- read.tree("Vascular_Plants_rooted.dated.tre")

library(readxl)
sp3 <- read_excel("Forest_data.xlsx", sheet = 6, range = "A2:A194", col_names = F)[[1]]

herb.list <- sapply(strsplit(sp3, " "), function(x) paste(x[1], x[2]))
herb.list <- sub(" ", "_", herb.list)

ha <- read_excel("Forest_data.xlsx", sheet = 6, range = "B2:CS194", col_names = F)
ha <- t(ha)
colnames(ha) <- herb.list
ha <- ha[, colSums(ha) > 0]
ha <- data.frame(ha)

source("addition-tools.r")

#_______________________________________________________

ha$Asclepiadaceae_Sp. <- NULL
ha$Caryophyllaceae_Sp. <- NULL
ha$Compositae_Sp. <- NULL
ha$Compositae_Sp..1 <- NULL
ha$Compositae_Sp..2 <- NULL
ha$Compositae_Sp..3 <- NULL
ha$Convolvulaceae_Sp. <- NULL
ha$Gramineae_Sp. <- NULL
ha$Gramineae_Sp..1 <- NULL
ha$Gramineae_Sp..2 <- NULL
ha$Gramineae_Sp..3 <- NULL
ha$Gramineae_Sp..4 <- NULL
ha$Gramineae_Sp..5 <- NULL
ha$Spp1.unclassified_NA <- NULL
ha$Spp2.unclassified_NA <- NULL
ha$Spp3.unclassified_NA <- NULL
ha$Umbelliferae_Sp. <- NULL

sp.list2 <- colnames(ha)

genera <- sapply(strsplit(colnames(ha), "_"), function(x) paste(x[1]))
g <- table(genera)

id2 <- which(sp.list2 %in% phylo$tip.label)
sp.list3 <- sp.list2[-id2]

#_______________________________________________________

grep("Bidens*", phylo$tip.label, value = T)

#Misprints
colnames(ha)[colnames(ha) == "Leibitzia_anandria"] <- "Leibnitzia_anandria"
colnames(ha)[colnames(ha) == "Impatiens_noli.tangere"] <- "Impatiens_noli"
colnames(ha)[colnames(ha) == "Lilium_lanciforlium"] <- "Lilium_lancifolium"
colnames(ha)[colnames(ha) == "Potenilla_supina"] <- "Potentilla_supina"
colnames(ha)[colnames(ha) == "Ploygonatum_odoratum"] <- "Polygonatum_odoratum"
colnames(ha)[colnames(ha) == "Tataxacum_mongolicum"] <- "Taraxacum_mongolicum"
colnames(ha)[colnames(ha) == "Oxyropis_coerulea"] <- "Oxytropis_caerulea"
colnames(ha)[colnames(ha) == "Patrinia_scabiosaefolia"] <- "Patrinia_scabiosifolia"

one_sp <- c()
for (ii in 1:length(sp.list3)){
  #  one_sp <- c(one_sp, grep(names(g[as.vector(g) == 1])[ii], sp.list3, value = T))
  one_sp <- c(one_sp, grep(gsub("\\_.*","", sp.list3)[ii], names(g[as.vector(g) == 1]), value = T))
}

#Monotypic genera
phylo$tip.label[phylo$tip.label == "Amphicarpaea_bracteata"] <- "Amphicarpaea_trisperma"
phylo$tip.label[phylo$tip.label == "Aquilegia_viscosa"] <- "Aquilegia_yabeana"
phylo$tip.label[phylo$tip.label == "Bidens_mitis"] <- "Bidens_parviflora"
phylo$tip.label[phylo$tip.label == "Parasenecio_hastiformis"] <- "Cacalia_hastata"
phylo$tip.label[phylo$tip.label == "Chenopodium_urbicum"] <- "Chenopodium_giaucum"
phylo$tip.label[phylo$tip.label == "Comarum_palustre"] <- "Comarum_gracillima"
phylo$tip.label[phylo$tip.label == "Corydalis_flavula"] <- "Corydalis_raddeana"
phylo$tip.label[phylo$tip.label == "Cyperus_strigosus"] <- "Cyperus_Sp."
phylo$tip.label[phylo$tip.label == "Chrysanthemum_indicum"] <- "Dendranthema_lavandulifolium"
phylo$tip.label[phylo$tip.label == "Dracocephalum_nutans"] <- "Dracocephalum_rupestre"
phylo$tip.label[phylo$tip.label == "Erigeron_pulchellus"] <- "Erigeron_acer"
phylo$tip.label[phylo$tip.label == "Erysimum_cheiri"] <- "Erysimum_bungei"
phylo$tip.label[phylo$tip.label == "Aster_hispidus"] <- "Heteropappus_hispidus"
phylo$tip.label[phylo$tip.label == "Lespedeza_juncea"] <- "Lespedeza_formosa"
phylo$tip.label[phylo$tip.label == "Lychnis_coronaria"] <- "Lychnis_fulgens"
phylo$tip.label[phylo$tip.label == "Lysimachia_candida"] <- "Lysimachia_pentapetala"
phylo$tip.label[phylo$tip.label == "Mentha_spicata"] <- "Mentha_haplocalyx"
phylo$tip.label[phylo$tip.label == "Orobanche_hederae"] <- "Orobanche_coerulescens"
phylo$tip.label[phylo$tip.label == "Pedicularis_spicata"] <- "Pedicularis_striata"
phylo$tip.label[phylo$tip.label == "Pseudostellaria_jamesiana"] <- "Pseudostellaria_davidii"
phylo$tip.label[phylo$tip.label == "Isodon_coetsa"] <- "Rabdosia_japonica"
phylo$tip.label[phylo$tip.label == "Rhaponticum_repens"] <- "Rhaponticum_uniflorum"
phylo$tip.label[phylo$tip.label == "Rheum_spiciforme"] <- "Rheum_franzenbachii"
phylo$tip.label[phylo$tip.label == "Sedum_rupestre"] <- "Sedum_aizoon"
phylo$tip.label[phylo$tip.label == "Stellaria_nemorum"] <- "Stellaria_dichotoma"
phylo$tip.label[phylo$tip.label == "Thesium_chinense"] <- "Thesium_refractum"
phylo$tip.label[phylo$tip.label == "Gentiana_lutea"] <- "Tripterospermum_chinense"
phylo$tip.label[phylo$tip.label == "Veronica_lanceolata"] <- "Veronica_linariifolia"

#Polytypic genera
phylo <- my.add(phylo, "Astragalus_sp.", tol = 1e-6)
phylo <- my.add(phylo, "Iris_Sp.", tol = 1e-6)
phylo <- my.add(phylo, "Rubus_Sp.", tol = 1e-6)
phylo <- my.add(phylo, "Taraxacum_Sp.", tol = 1e-6)
phylo <- my.add(phylo, "Vicia_Sp.", tol = 1e-6)
phylo <- my.add(phylo, "Viola_Sp.", tol = 1e-6)

#Closest species
phylo$tip.label[phylo$tip.label == "Aconitum_sinomontanum"] <- "Aconitum_barbatum" #doi.org/10.1371/journal.pone.0171038

phylo$tip.label[phylo$tip.label == "Adenophora_triphylla"] <- "Adenophora_divaricata" #DOI: 10.1600/036364408783887465
phylo$tip.label[phylo$tip.label == "Adenophora_stricta"] <- "Adenophora_polyantha"
phylo$tip.label[phylo$tip.label == "Adenophora_verticillata"] <- "Adenophora_wawreana"

phylo$tip.label[phylo$tip.label == "Artemisia_laciniata"] <- "Artemisia_tanacetifolia"  #Smith, 2011

phylo$tip.label[phylo$tip.label == "Clematis_villosa"] <- "Clematis_ochotensis" #Smith, 2011

phylo$tip.label[phylo$tip.label == "Ixeris_repens"] <- "Ixeris_polycephala" #doi.org/10.1508/cytologia.FujiiJubilaei.188

phylo$tip.label[phylo$tip.label == "Patrinia_heterophylla"] <- "Patrinia_scabra" #Smith, 2011

phylo$tip.label[phylo$tip.label == "Polygonatum cirrhifolium"] <- "Polygonatum_sibiricum" #Smith, 2011

phylo$tip.label[phylo$tip.label == "Polygonum_weyrichii"] <- "Polygonum_bistorta" #Smith, 2011

phylo$tip.label[phylo$tip.label == "Prenanthes_autumnalis"] <- "Prenanthes_tatarinowii"

phylo$tip.label[phylo$tip.label == "Vicia_unijuga"] <- "Vicia_pseudo-orobus" #doi.org/10.1186/1471-2148-12-250

phylo$tip.label[phylo$tip.label == "Viola_lactiflora"] <- "Viola_yezoensis"

#

ha$Prenanthes_tatarinowii <- ha$Prenanthes_macrophylla + ha$Prenanthes_tatarinowii #the same species
ha$Prenanthes_macrophylla <- NULL

ha$Bupleurum_chinense <- ha$Bupleurum_chinense + ha$Bupleurum_chinense.1 #Subspecies
ha$Bupleurum_chinense.1 <- NULL

ha$Thalictrum_minus <- ha$Thalictrum_minus + ha$Thalictrum_var.
ha$Thalictrum_var. <- NULL

#The rest of species

sp.list2 <- colnames(ha)
id2 <- which(sp.list2 %in% phylo$tip.label)
sp.list3 <- sp.list2[-id2]

for (ii in 1:length(sp.list3)){
  phylo <- my.add(phylo, sp.list3[ii], tol = 1e-6)
}

colnames(ha)[colnames(ha) == "Denfranthema_chanetii"] <- "Dendranthema_chanetii"
phylo <- my.add(phylo, "Dendranthema_chanetii", tol = 1e-6)

#_______________________________________________________

sp.list2 <- colnames(ha)
id <- which(phylo$tip.label %in% sp.list2)
tree <- drop.tip(phylo, phylo$tip.label[-id])
plot(tree, cex = 0.7)
axisPhylo()

rm(g, genera, herb.list, id, id2, ii, one_sp, sp.list2, sp.list3, sp3)
