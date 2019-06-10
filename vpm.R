devtools::install_github("jinyizju/V.PhyloMaker")

library(V.PhyloMaker)


#таблицы с представленностями, все дела

library(readxl)
sp1 <- read_excel("raw.data/Forest_data.xlsx", sheet = 4, range = "A2:A21", col_names = F)[[1]]
sp2 <- read_excel("raw.data/Forest_data.xlsx", sheet = 5, range = "A2:A42", col_names = F)[[1]]

tree.list <- sapply(strsplit(sp1, " "), function(x) paste(x[1], x[2]))
tree.list <- sub(" ", "_", tree.list)

shrub.list <- sapply(strsplit(sp2, " "), function(x) paste(x[1], x[2]))
shrub.list <- sub(" ", "_", shrub.list)

ta <- read_excel("raw.data/Forest_data.xlsx", sheet = 4, range = "B2:CS21", col_names = F)
ta <- t(ta)
colnames(ta) <- tree.list
ta <- ta[, colSums(ta) > 0]
ta <- data.frame(ta)

sa <- read_excel("raw.data/Forest_data.xlsx", sheet = 5, range = "B2:CS42", col_names = F)
sa <- t(sa)
colnames(sa) <- shrub.list
sa <- sa[, colSums(sa) > 0]
sa <- data.frame(sa)

#_________________________________________________________________

#Функция требует список видов в форме таблицы, состоящей из 5 столбцов:
#вид, род, семейство, ближайший вид, ближайший род (последние два могут быть пустыми, 
#они используются при добавлении видов и родов, отсутствующих в филогении)
sp.list <- read.csv2("species list.csv")

#По умолчанию функция использует древо GBOTB.extended, но может быть использовано любое другое.
#В этом случае необходимо при помощи функции build.nodes.1/build.nodes.2, которая извлекает информацию
#о расстоянии, соединяющих узлы в дереве. Я так и не поняла разницу между ними, но 
#функции по-разному работают с не монофилетическими родами. 

#nodes запрашивает результат функции build.nodes, для дерева GBOTB.extended они уже присутствуют в пакете - 
#nodes.info.1 и nodes.info.2 соответственно.

#scenarios определяет, каким образом будут добавляться новые вершины в дерево
#S1 - новая вершина привязывается к основанию рода/семейства (как политомия, судя по всему)
#S2 - привязывается к случайному узлу на уровне рода/семейства
#S3 - привязывается к середине ветви рода/семейства (как я поняла)
tree.a <- phylo.maker(sp.list = sp.list, tree = GBOTB.extended, nodes = nodes.info.1, scenarios="S1")

#plot(tree.a$scenario.1, cex = 0.7)
#axisPhylo()

#_________________________________________________________________

#это все можно не выполнять

wood.tree <- tree.a$scenario.1

library(picante)
wood.mat <- cophenetic(tree.a$scenario.1)

t.ses.tl <- ses.mpd(ta, wood.mat, null.model = "taxa.labels", runs = 999, iterations = 1000)
t.ses.r  <- ses.mpd(ta, wood.mat, null.model = "richness", runs = 999, iterations = 1000)
t.ses.f  <- ses.mpd(ta, wood.mat, null.model = "frequency", runs = 999, iterations = 1000)
t.ses.sp <- ses.mpd(ta, wood.mat, null.model = "sample.pool", runs = 999, iterations = 1000)
t.ses.pp <- ses.mpd(ta, wood.mat, null.model = "phylogeny.pool", runs = 999, iterations = 1000)
t.ses.is <- ses.mpd(ta, wood.mat, null.model = "independentswap", runs = 999, iterations = 1000)
t.ses.ts <- ses.mpd(ta, wood.mat, null.model = "trialswap", runs = 999, iterations = 1000)

t.ses.a.tl <- ses.mpd(ta, wood.mat, null.model = "taxa.labels", runs = 999, iterations = 1000, abundance.weighted = T)
t.ses.a.r  <- ses.mpd(ta, wood.mat, null.model = "richness", runs = 999, iterations = 1000, abundance.weighted = T)
t.ses.a.f  <- ses.mpd(ta, wood.mat, null.model = "frequency", runs = 999, iterations = 1000, abundance.weighted = T)
t.ses.a.sp <- ses.mpd(ta, wood.mat, null.model = "sample.pool", runs = 999, iterations = 1000, abundance.weighted = T)
t.ses.a.pp <- ses.mpd(ta, wood.mat, null.model = "phylogeny.pool", runs = 999, iterations = 1000, abundance.weighted = T)
t.ses.a.is <- ses.mpd(ta, wood.mat, null.model = "independentswap", runs = 999, iterations = 1000, abundance.weighted = T)
t.ses.a.ts <- ses.mpd(ta, wood.mat, null.model = "trialswap", runs = 999, iterations = 1000, abundance.weighted = T)

hgt <- read.table("raw.data/hgt.txt")[[1]]

t.ses.z <- data.frame(hgt = hgt, 
                      taxa.labels = t.ses.tl$mpd.obs.z,
                      richness = t.ses.r$mpd.obs.z,
                      frequency = t.ses.f$mpd.obs.z,
                      sample.pool = t.ses.sp$mpd.obs.z,
                      phylogeny.pool = t.ses.pp$mpd.obs.z,
                      independentswap = t.ses.is$mpd.obs.z,
                      trialswap = t.ses.ts$mpd.obs.z,
                      taxa.labels.a = t.ses.a.tl$mpd.obs.z, 
                      richness.a = t.ses.a.r$mpd.obs.z,
                      frequency.a = t.ses.a.f$mpd.obs.z,
                      sample.pool.a = t.ses.a.sp$mpd.obs.z,
                      phylogeny.pool.a = t.ses.a.pp$mpd.obs.z,
                      independentswap.a = t.ses.a.is$mpd.obs.z,
                      trialswap.a = t.ses.a.ts$mpd.obs.z)

t.ses.p <- data.frame(hgt = hgt, 
                      taxa.labels = t.ses.tl$mpd.obs.p,
                      richness = t.ses.r$mpd.obs.p,
                      frequency = t.ses.f$mpd.obs.p,
                      sample.pool = t.ses.sp$mpd.obs.p,
                      phylogeny.pool = t.ses.pp$mpd.obs.p,
                      independentswap = t.ses.is$mpd.obs.p,
                      trialswap = t.ses.ts$mpd.obs.p,
                      taxa.labels.a = t.ses.a.tl$mpd.obs.p, 
                      richness.a = t.ses.a.r$mpd.obs.p,
                      frequency.a = t.ses.a.f$mpd.obs.p,
                      sample.pool.a = t.ses.a.sp$mpd.obs.p,
                      phylogeny.pool.a = t.ses.a.pp$mpd.obs.p,
                      independentswap.a = t.ses.a.is$mpd.obs.p,
                      trialswap.a = t.ses.a.ts$mpd.obs.p)

save(t.ses.z, t.ses.p, file = "vpm-phylo-ses.rda")

#_________________________________________________________________

#а вот это можно выполнить

load("vpm-phylo-ses.rda")
hgt <- read.table("raw.data/hgt.txt")[[1]]

pdf("figures/tree-regs2.pdf", 7, 22)

op <- par(mfcol = c(7,2))

for (i in 2:ncol(t.ses.z)){
  NRI <- -t.ses.z[, i]
  plot(hgt, NRI,  pch = 21, bg = "tomato", xlab = "altitude")
  
  sp <- cor.test(hgt, NRI, method = "spearman", use = "complete")
  
  sig <- t.ses.p[, i] < 0.05 | t.ses.p[, i] > 0.95
  points(hgt[sig], NRI[sig],  pch = 21, bg = "darkred")
  
  model <- lm(NRI ~ hgt)
  if (sp$p.value < 0.05)
  {
    abline(model)
  }
  legend("topright", legend = round(coef(model)[2]*100, 4), bty = "n")
  
  text("topleft", labels = "txt")
  legend("topleft", legend = bquote(rho == .({round(sp$estimate, 3)}) ~ "; " ~ p == .({round(sp$p.value, 3)})), bty = "n")
  
  title(main = paste0(names(t.ses.z)[i]))
}

par(op)

dev.off()

#Графики получились очень похожими на те, что мы делали по старой схеме.
