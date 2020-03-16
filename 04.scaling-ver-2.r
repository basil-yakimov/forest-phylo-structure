load("clean.data/tree-abund.rda")
load("clean.data/shrub-abund.rda")
load("clean.data/herb-abund.rda")

load("clean.data/wood-phylo.rda")
load("clean.data/herb-phylo.rda")


source("R/scaling-ses.r")



t.dist <- cophenetic(drop.tip(wood.tree, wood.tree$tip.label[!(wood.tree$tip.label %in% colnames(ta))]))
t.dist <- t.dist[names(ta), names(ta)]

t.grain <- scaling.grain(ta, t.dist, iter = 1000)
t.grain.a <- scaling.grain(ta, t.dist, ab = T, iter = 1000)

t.pool <- scaling.pool(ta, t.dist, iter = 1000)
t.pool.a <- scaling.pool(ta, t.dist, ab = T, iter = 1000)

t.coherent <- scaling.coherent(ta, t.dist, iter = 1000)
t.coherent.a <- scaling.coherent(ta, t.dist, ab = T, iter = 1000)



s.dist <- cophenetic(drop.tip(wood.tree, wood.tree$tip.label[!(wood.tree$tip.label %in% colnames(sa))]))
s.dist <- s.dist[names(sa), names(sa)]

s.grain <- scaling.grain(sa, s.dist, iter = 1000)
s.grain.a <- scaling.grain(sa, s.dist, ab = T, iter = 1000)

s.pool <- scaling.pool(sa, s.dist, iter = 1000)
s.pool.a <- scaling.pool(sa, s.dist, ab = T, iter = 1000)

s.coherent <- scaling.coherent(sa, s.dist, iter = 1000)
s.coherent.a <- scaling.coherent(sa, s.dist, ab = T, iter = 1000)



h.dist <- cophenetic(drop.tip(herb.tree, herb.tree$tip.label[!(herb.tree$tip.label %in% colnames(ha))]))
h.dist <- h.dist[names(ha), names(ha)]

h.grain <- scaling.grain(ha, h.dist, iter = 1000)
h.grain.a <- scaling.grain(ha, h.dist, ab = T, iter = 1000)

h.pool <- scaling.pool(ha, h.dist, iter = 1000)
h.pool.a <- scaling.pool(ha, h.dist, ab = T, iter = 1000)

h.coherent <- scaling.coherent(ha, h.dist, iter = 1000)
h.coherent.a <- scaling.coherent(ha, h.dist, ab = T, iter = 1000)


save(t.grain, t.grain.a,
     s.grain, s.grain.a, 
     h.grain, h.grain.a, file = "clean.data/grain-scaling.rda")

save(t.pool, t.pool.a,
     s.pool, s.pool.a,
     h.pool, h.pool.a, file = "clean.data/pool-scaling.rda")

save(t.coherent, t.coherent.a,
     s.coherent, s.coherent.a,
     h.coherent, h.coherent.a, file = "clean.data/coherent-scaling.rda")


#system('shutdown -s')

#_______________________________________________________________________________________________

#functional

load("clean.data/func-abund.rda")

load("clean.data/tree-traits.rda")
load("clean.data/shrub-traits.rda")
load("clean.data/herb-traits.rda")


source("R/scaling-ses.r")



#t.dist <- t.dist[names(ta), names(ta)]
t.dist <- t.dist[order(rownames(t.dist)), order(colnames(t.dist))]

t.grain <- scaling.grain(ta, t.dist, iter = 1000)
t.grain.a <- scaling.grain(ta, t.dist, ab = T, iter = 1000)

t.pool <- scaling.pool(ta, t.dist, iter = 1000)
t.pool.a <- scaling.pool(ta, t.dist, ab = T, iter = 1000)

t.coherent <- scaling.coherent(ta, t.dist, iter = 1000)
t.coherent.a <- scaling.coherent(ta, t.dist, ab = T, iter = 1000)



#s.dist <- s.dist[names(sa), names(sa)]
s.dist <- s.dist[order(rownames(s.dist)), order(colnames(s.dist))]

s.grain <- scaling.grain(sa, s.dist, iter = 1000)
s.grain.a <- scaling.grain(sa, s.dist, ab = T, iter = 1000)

s.pool <- scaling.pool(sa, s.dist, iter = 1000)
s.pool.a <- scaling.pool(sa, s.dist, ab = T, iter = 1000)

s.coherent <- scaling.coherent(sa, s.dist, iter = 1000)
s.coherent.a <- scaling.coherent(sa, s.dist, ab = T, iter = 1000)



#h.dist <- h.dist[names(ha), names(ha)]
h.dist <- h.dist[order(rownames(h.dist)), order(colnames(h.dist))]

h.grain <- scaling.grain(ha, h.dist, iter = 1000)
h.grain.a <- scaling.grain(ha, h.dist, ab = T, iter = 1000)

h.pool <- scaling.pool(ha, h.dist, iter = 1000)
h.pool.a <- scaling.pool(ha, h.dist, ab = T, iter = 1000)

h.coherent <- scaling.coherent(ha, h.dist, iter = 1000)
h.coherent.a <- scaling.coherent(ha, h.dist, ab = T, iter = 1000)

save(t.grain, t.grain.a,
     s.grain, s.grain.a, 
     h.grain, h.grain.a, file = "clean.data/func-grain-scaling.rda")

save(t.pool, t.pool.a,
     s.pool, s.pool.a,
     #h.pool, h.pool.a,
     file = "clean.data/func-pool-scaling.rda")

save(t.coherent, t.coherent.a,
     s.coherent, s.coherent.a,
     #h.coherent, h.coherent.a,
     file = "clean.data/func-coherent-scaling.rda")


#system('shutdown -s')
