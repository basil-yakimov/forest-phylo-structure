load("clean.data/grain-scaling.rda")
load("clean.data/pool-scaling.rda")
load("clean.data/coherent-scaling.rda")

source("R/scaling-ses.r")


png("figures/tree-phylo-grain-scaling.png", width = 700, height = 500)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot.scaling(t.grain, col = "#ffbcb0")
dev.off()

png("figures/shrub-phylo-grain-scaling.png", width = 700, height = 500)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot.scaling(s.grain, col = "#b6eaff")
dev.off()

png("figures/herb-phylo-grain-scaling.png", width = 700, height = 500)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot.scaling(h.grain, col = "#b3ffb3")
dev.off()



png("figures/tree-phylo-grain-scaling-a.png", width = 700, height = 500)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot.scaling(t.grain.a, col = "#ffbcb0")
dev.off()

png("figures/shrub-phylo-grain-scaling-a.png", width = 700, height = 500)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot.scaling(s.grain.a, col = "#b6eaff")
dev.off()

png("figures/herb-phylo-grain-scaling-a.png", width = 700, height = 500)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot.scaling(h.grain.a, col = "#b3ffb3")
dev.off()

#_________________________________________________________________________________________________


png("figures/tree-phylo-pool-scaling.png", width = 700, height = 500)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot.scaling(t.pool, col = "#ffbcb0")
dev.off()

png("figures/shrub-phylo-pool-scaling.png", width = 700, height = 500)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot.scaling(s.pool, col = "#b6eaff")
dev.off()

png("figures/herb-phylo-pool-scaling.png", width = 700, height = 500)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot.scaling(h.pool, col = "#b3ffb3")
dev.off()



png("figures/tree-phylo-pool-scaling-a.png", width = 700, height = 500)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot.scaling(t.pool.a, col = "#ffbcb0")
dev.off()

png("figures/shrub-phylo-pool-scaling-a.png", width = 700, height = 500)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot.scaling(s.pool.a, col = "#b6eaff")
dev.off()

png("figures/herb-phylo-pool-scaling-a.png", width = 700, height = 500)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot.scaling(h.pool.a, col = "#b3ffb3")
dev.off()

#_________________________________________________________________________________________________


png("figures/tree-phylo-coherent-scaling.png", width = 700, height = 500)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot.scaling(t.coherent, col = "#ffbcb0")
dev.off()

png("figures/shrub-phylo-coherent-scaling.png", width = 700, height = 500)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot.scaling(s.coherent, col = "#b6eaff")
dev.off()

png("figures/herb-phylo-coherent-scaling.png", width = 700, height = 500)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot.scaling(h.coherent, col = "#b3ffb3")
dev.off()



png("figures/tree-phylo-coherent-scaling-a.png", width = 700, height = 500)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot.scaling(t.coherent.a, col = "#ffbcb0")
dev.off()

png("figures/shrub-phylo-coherent-scaling-a.png", width = 700, height = 500)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot.scaling(s.coherent.a, col = "#b6eaff")
dev.off()

png("figures/herb-phylo-coherent-scaling-a.png", width = 700, height = 500)
par(mar = c(4, 4, 0.5, 0.5), cex = 2)
plot.scaling(h.coherent.a, col = "#b3ffb3")
dev.off()
