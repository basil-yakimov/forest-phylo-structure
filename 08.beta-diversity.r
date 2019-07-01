load("clean.data/tree-abund.rda")
load("clean.data/shrub-abund.rda")
load("clean.data/herb-abund.rda")
load("clean.data/wood-phylo.rda")
load("clean.data/herb-phylo.rda")

hgt <- read.table("raw.data/hgt.txt")[[1]]

#1110-1180,1190-1255, 1260-1338, 1344-1408, 1414-1496, 1504-1552, 1558-1616, 1622-1682, and 1687-1734
subs <- cut(hgt, breaks = c(0, 1185, 1257, 1340, 1410, 1500, 1554, 1620, 1685, 3000),
            labels = paste0("sub", 1:9))

table(subs)

#________________________________________

source("R/facet-tools.r")

# tree layer taxonomic and phylogenetic beta-diversity

res0 <- tapply(1:96, subs, function(x) t.decomp(ta[x, ], 0), simplify = T)
res0 <- data.frame(t(simplify2array(res0)))

res1 <- tapply(1:96, subs, function(x) t.decomp(ta[x, ], 1), simplify = T)
res1 <- data.frame(t(simplify2array(res1)))

res2 <- tapply(1:96, subs, function(x) t.decomp(ta[x, ], 2), simplify = T)
res2 <- data.frame(t(simplify2array(res2)))

tb <- data.frame(c0 = res0$c, c1 = res1$c, c2 = res2$c,
                  u0 = res0$u, u1 = res1$u, u2 = res2$u)


res0 <- tapply(1:96, subs, function(x) p.decomp(ta[x, ], wood.tree, 0), simplify = T)
res0 <- data.frame(t(simplify2array(res0)))

res1 <- tapply(1:96, subs, function(x) p.decomp(ta[x, ], wood.tree, 1), simplify = T)
res1 <- data.frame(t(simplify2array(res1)))

res2 <- tapply(1:96, subs, function(x) p.decomp(ta[x, ], wood.tree, 2), simplify = T)
res2 <- data.frame(t(simplify2array(res2)))

tb <- cbind(tb, data.frame(pc0 = res0$c, pc1 = res1$c, pc2 = res2$c,
                            pu0 = res0$u, pu1 = res1$u, pu2 = res2$u))

# shrub layer taxonomic and phylogenetic beta-diversity

res0 <- tapply(1:96, subs, function(x) t.decomp(sa[x, ], 0), simplify = T)
res0 <- data.frame(t(simplify2array(res0)))

res1 <- tapply(1:96, subs, function(x) t.decomp(sa[x, ], 1), simplify = T)
res1 <- data.frame(t(simplify2array(res1)))

res2 <- tapply(1:96, subs, function(x) t.decomp(sa[x, ], 2), simplify = T)
res2 <- data.frame(t(simplify2array(res2)))

sb <- data.frame(c0 = res0$c, c1 = res1$c, c2 = res2$c,
                 u0 = res0$u, u1 = res1$u, u2 = res2$u)


res0 <- tapply(1:96, subs, function(x) p.decomp(sa[x, ], wood.tree, 0), simplify = T)
res0 <- data.frame(t(simplify2array(res0)))

res1 <- tapply(1:96, subs, function(x) p.decomp(sa[x, ], wood.tree, 1), simplify = T)
res1 <- data.frame(t(simplify2array(res1)))

res2 <- tapply(1:96, subs, function(x) p.decomp(sa[x, ], wood.tree, 2), simplify = T)
res2 <- data.frame(t(simplify2array(res2)))

sb <- cbind(sb, data.frame(pc0 = res0$c, pc1 = res1$c, pc2 = res2$c,
                           pu0 = res0$u, pu1 = res1$u, pu2 = res2$u))

# herb layer taxonomic and phylogenetic beta-diversity

res0 <- tapply(1:96, subs, function(x) t.decomp(ha[x, ], 0), simplify = T)
res0 <- data.frame(t(simplify2array(res0)))

res1 <- tapply(1:96, subs, function(x) t.decomp(ha[x, ], 1), simplify = T)
res1 <- data.frame(t(simplify2array(res1)))

res2 <- tapply(1:96, subs, function(x) t.decomp(ha[x, ], 2), simplify = T)
res2 <- data.frame(t(simplify2array(res2)))

hb <- data.frame(c0 = res0$c, c1 = res1$c, c2 = res2$c,
                 u0 = res0$u, u1 = res1$u, u2 = res2$u)


res0 <- tapply(1:96, subs, function(x) p.decomp(ha[x, ], herb.tree, 0), simplify = T)
res0 <- data.frame(t(simplify2array(res0)))

res1 <- tapply(1:96, subs, function(x) p.decomp(ha[x, ], herb.tree, 1), simplify = T)
res1 <- data.frame(t(simplify2array(res1)))

res2 <- tapply(1:96, subs, function(x) p.decomp(ha[x, ], herb.tree, 2), simplify = T)
res2 <- data.frame(t(simplify2array(res2)))

hb <- cbind(hb, data.frame(pc0 = res0$c, pc1 = res1$c, pc2 = res2$c,
                           pu0 = res0$u, pu1 = res1$u, pu2 = res2$u))


#-----------------------------------------#

source("R/plot.ses.r")
mhgt <- tapply(hgt, subs, mean)


op <- par(mfrow = c(3, 6), mar = c(4, 4, 0.5, 0.5))
plot.ses(tb$c0, mhgt, col = "forestgreen", lab = expression(C[0]))
plot.ses(tb$pc0, mhgt, col = "forestgreen", lab = expression(PC[0]))
plot.ses(sb$c0, mhgt, col = "skyblue", lab = expression(C[0]))
plot.ses(sb$pc0, mhgt, col = "skyblue", lab = expression(PC[0]))
plot.ses(tb$c0, mhgt, col = "tomato", lab = expression(C[0]))
plot.ses(tb$pc0, mhgt, col = "tomato", lab = expression(PC[0]))

plot.ses(tb$c1, mhgt, col = "forestgreen", lab = expression(C[1]))
plot.ses(tb$pc1, mhgt, col = "forestgreen", lab = expression(PC[1]))
plot.ses(sb$c1, mhgt, col = "skyblue", lab = expression(C[1]))
plot.ses(sb$pc1, mhgt, col = "skyblue", lab = expression(PC[1]))
plot.ses(tb$c1, mhgt, col = "tomato", lab = expression(C[1]))
plot.ses(tb$pc1, mhgt, col = "tomato", lab = expression(PC[1]))

plot.ses(tb$c2, mhgt, col = "forestgreen", lab = expression(C[2]))
plot.ses(tb$pc2, mhgt, col = "forestgreen", lab = expression(PC[2]))
plot.ses(sb$c2, mhgt, col = "skyblue", lab = expression(C[2]))
plot.ses(sb$pc2, mhgt, col = "skyblue", lab = expression(PC[2]))
plot.ses(tb$c2, mhgt, col = "tomato", lab = expression(C[2]))
plot.ses(tb$pc2, mhgt, col = "tomato", lab = expression(PC[2]))
par(op)

op <- par(mfrow = c(3, 6), mar = c(4, 4, 0.5, 0.5))
plot.ses(tb$u0, mhgt, col = "forestgreen", lab = expression(U[0]))
plot.ses(tb$pu0, mhgt, col = "forestgreen", lab = expression(PU[0]))
plot.ses(sb$u0, mhgt, col = "skyblue", lab = expression(U[0]))
plot.ses(sb$pu0, mhgt, col = "skyblue", lab = expression(PU[0]))
plot.ses(tb$u0, mhgt, col = "tomato", lab = expression(U[0]))
plot.ses(tb$pu0, mhgt, col = "tomato", lab = expression(PU[0]))

plot.ses(tb$u1, mhgt, col = "forestgreen", lab = expression(U[1]))
plot.ses(tb$pu1, mhgt, col = "forestgreen", lab = expression(PU[1]))
plot.ses(sb$u1, mhgt, col = "skyblue", lab = expression(U[1]))
plot.ses(sb$pu1, mhgt, col = "skyblue", lab = expression(PU[1]))
plot.ses(tb$u1, mhgt, col = "tomato", lab = expression(U[1]))
plot.ses(tb$pu1, mhgt, col = "tomato", lab = expression(PU[1]))

plot.ses(tb$u2, mhgt, col = "forestgreen", lab = expression(U[2]))
plot.ses(tb$pu2, mhgt, col = "forestgreen", lab = expression(PU[2]))
plot.ses(sb$u2, mhgt, col = "skyblue", lab = expression(U[2]))
plot.ses(sb$pu2, mhgt, col = "skyblue", lab = expression(PU[2]))
plot.ses(tb$u2, mhgt, col = "tomato", lab = expression(U[2]))
plot.ses(tb$pu2, mhgt, col = "tomato", lab = expression(PU[2]))
par(op)

