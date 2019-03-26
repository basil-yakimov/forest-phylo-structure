load("clean.data/scaling.rda")
load("clean.data/scaling-ses.rda")

#_____________________________________________________________________


png("figures/regs-sc.png", 4000, 24000, pointsize = 75)

op <- par(mfrow = c(20,3))

for (ii in 1:20)
{
    hgt <- hgt_sc[[ii]]
    
    NRI <- -ta_sc_ses[[ii]][, "z.is"]
    plot(hgt, NRI,  pch = 21, bg = "tomato", xlab = "altitude")
    model <- lm(NRI ~ hgt)
    if (anova(model)[1, 5] < 0.05) abline(model)
    sp <- cor.test(hgt, NRI, method = "spearman", use = "complete")
    title(bquote("Tree " ~ .(ii) ~ ": " ~ rho == .({round(sp$estimate, 3)}) ~ "; " ~ p == .({round(sp$p.value, 3)}) ~ "; " ~ R^2 == .({round(cor(hgt, NRI, use = "complete")^2, 3)})))

    NRI <- -sa_sc_ses[[ii]][, "z.is"]
    plot(hgt, NRI,  pch = 21, bg = "tomato", xlab = "altitude")
    model <- lm(NRI ~ hgt)
    if (anova(model)[1, 5] < 0.05) abline(model)
    sp <- cor.test(hgt, NRI, method = "spearman", use = "complete")
    title(bquote("Shrub " ~ .(ii) ~ ": " ~ rho == .({round(sp$estimate, 3)}) ~ "; " ~ p == .({round(sp$p.value, 3)}) ~ "; " ~ R^2 == .({round(cor(hgt, NRI, use = "complete")^2, 3)})))
    
    NRI <- -ha_sc_ses[[ii]][, "z.is"]
    plot(hgt, NRI,  pch = 21, bg = "tomato", xlab = "altitude")
    model <- lm(NRI ~ hgt)
    if (anova(model)[1, 5] < 0.05) abline(model)
    sp <- cor.test(hgt, NRI, method = "spearman", use = "complete")
    title(bquote("Herb " ~ .(ii) ~ ": " ~ rho == .({round(sp$estimate, 3)}) ~ "; " ~ p == .({round(sp$p.value, 3)}) ~ "; " ~ R^2 == .({round(cor(hgt, NRI, use = "complete")^2, 3)})))
}

dev.off()

png("figures/regs-sc-a.png", 4000, 24000, pointsize = 75)

op <- par(mfrow = c(20,3))

for (ii in 1:20)
{
  hgt <- hgt_sc[[ii]]
  
  NRI <- -ta_sc_ses[[ii]][, "z.a.is"]
  plot(hgt, NRI,  pch = 21, bg = "tomato", xlab = "altitude")
  model <- lm(NRI ~ hgt)
  if (anova(model)[1, 5] < 0.05) abline(model)
  sp <- cor.test(hgt, NRI, method = "spearman", use = "complete")
  title(bquote("Tree " ~ .(ii) ~ ": " ~ rho == .({round(sp$estimate, 3)}) ~ "; " ~ p == .({round(sp$p.value, 3)}) ~ "; " ~ R^2 == .({round(cor(hgt, NRI, use = "complete")^2, 3)})))
  
  NRI <- -sa_sc_ses[[ii]][, "z.a.is"]
  plot(hgt, NRI,  pch = 21, bg = "tomato", xlab = "altitude")
  model <- lm(NRI ~ hgt)
  if (anova(model)[1, 5] < 0.05) abline(model)
  sp <- cor.test(hgt, NRI, method = "spearman", use = "complete")
  title(bquote("Shrub " ~ .(ii) ~ ": " ~ rho == .({round(sp$estimate, 3)}) ~ "; " ~ p == .({round(sp$p.value, 3)}) ~ "; " ~ R^2 == .({round(cor(hgt, NRI, use = "complete")^2, 3)})))
  
  NRI <- -ha_sc_ses[[ii]][, "z.a.is"]
  plot(hgt, NRI,  pch = 21, bg = "tomato", xlab = "altitude")
  model <- lm(NRI ~ hgt)
  if (anova(model)[1, 5] < 0.05) abline(model)
  sp <- cor.test(hgt, NRI, method = "spearman", use = "complete")
  title(bquote("Herb " ~ .(ii) ~ ": " ~ rho == .({round(sp$estimate, 3)}) ~ "; " ~ p == .({round(sp$p.value, 3)}) ~ "; " ~ R^2 == .({round(cor(hgt, NRI, use = "complete")^2, 3)})))
}

dev.off()


#_____________________________________________________________________

sc <- 1:20

t.res <- sapply(sc, function(x) {
  NRI <- -ta_sc_ses[[x]][, "z.is"]
  hgt <- hgt_sc[[x]]
  return(c(cor(NRI, hgt, use = "complete")^2, cor(NRI, hgt, method = "spearman", use = "complete")))
})

s.res <- sapply(sc, function(x) {
  NRI <- -sa_sc_ses[[x]][, "z.is"]
  hgt <- hgt_sc[[x]]
  return(c(cor(NRI, hgt, use = "complete")^2, cor(NRI, hgt, method = "spearman", use = "complete")))
})

h.res <- sapply(sc, function(x) {
  NRI <- -ha_sc_ses[[x]][, "z.is"]
  hgt <- hgt_sc[[x]]
  return(c(cor(NRI, hgt, use = "complete")^2, cor(NRI, hgt, method = "spearman", use = "complete")))
})


png("figures/coef-sc.png", 3000, 4000, pointsize = 75)

op <- par(mfrow = c(3,2))

plot(sc, t.res[2, ], pch = 21, bg = "tomato", ylab = expression(rho), main = "Tree")
plot(sc, t.res[1, ], pch = 21, bg = "tomato", ylab = expression(r^2), main = "Tree")

plot(sc, s.res[2, ], pch = 21, bg = "skyblue", ylab = expression(rho), main = "Shrub")
plot(sc, s.res[1, ], pch = 21, bg = "skyblue", ylab = expression(r^2), main = "Shrub")

plot(sc, h.res[2, ], pch = 21, bg = "forestgreen", ylab = expression(rho), main = "Herb")
plot(sc, h.res[1, ], pch = 21, bg = "forestgreen", ylab = expression(r^2), main = "Herb")

dev.off()



t.res <- sapply(sc, function(x) {
  NRI <- -ta_sc_ses[[x]][, "z.a.is"]
  hgt <- hgt_sc[[x]]
  return(c(cor(NRI, hgt, use = "complete")^2, cor(NRI, hgt, method = "spearman", use = "complete")))
})

s.res <- sapply(sc, function(x) {
  NRI <- -sa_sc_ses[[x]][, "z.a.is"]
  hgt <- hgt_sc[[x]]
  return(c(cor(NRI, hgt, use = "complete")^2, cor(NRI, hgt, method = "spearman", use = "complete")))
})

h.res <- sapply(sc, function(x) {
  NRI <- -ha_sc_ses[[x]][, "z.a.is"]
  hgt <- hgt_sc[[x]]
  return(c(cor(NRI, hgt, use = "complete")^2, cor(NRI, hgt, method = "spearman", use = "complete")))
})


png("figures/coef-a-sc.png", 3000, 4000, pointsize = 75)

op <- par(mfrow = c(3,2))

plot(sc, t.res[2, ], pch = 21, bg = "tomato", ylab = expression(rho), main = "Tree")
plot(sc, t.res[1, ], pch = 21, bg = "tomato", ylab = expression(r^2), main = "Tree")

plot(sc, s.res[2, ], pch = 21, bg = "skyblue", ylab = expression(rho), main = "Shrub")
plot(sc, s.res[1, ], pch = 21, bg = "skyblue", ylab = expression(r^2), main = "Shrub")

plot(sc, h.res[2, ], pch = 21, bg = "forestgreen", ylab = expression(rho), main = "Herb")
plot(sc, h.res[1, ], pch = 21, bg = "forestgreen", ylab = expression(r^2), main = "Herb")

dev.off()

#_____________________________________________________________________

png("figures/mean-nri-sc.png", 4000, 4000, pointsize = 75)

op <- par(mfrow = c(4,3), mar = c(3, 3, 2, 0.5))

meanNRI <- sapply(ta_sc_ses, function(x) -mean(x[, "z.tl"], na.rm = T))
plot(sc, meanNRI, pch = 21, bg = "tomato", main = "Tree: taxa labels")
model <- lm(meanNRI ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model)

meanNRI <- sapply(sa_sc_ses, function(x) -mean(x[, "z.tl"], na.rm = T))
plot(sc, meanNRI, pch = 21, bg = "skyblue", main = "Shrub: taxa labels")
model <- lm(meanNRI ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model)

meanNRI <- sapply(ha_sc_ses, function(x) -mean(x[, "z.tl"], na.rm = T))
plot(sc, meanNRI, pch = 21, bg = "forestgreen", main = "Herb: taxa labels")
model <- lm(meanNRI ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model)

#---#

meanNRI <- sapply(ta_sc_ses, function(x) -mean(x[, "z.r"], na.rm = T))
plot(sc, meanNRI, pch = 21, bg = "tomato", main = "Tree: richness")
model <- lm(meanNRI ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model)

meanNRI <- sapply(sa_sc_ses, function(x) -mean(x[, "z.r"], na.rm = T))
plot(sc, meanNRI, pch = 21, bg = "skyblue", main = "Shrub: richness")
model <- lm(meanNRI ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model)

meanNRI <- sapply(ha_sc_ses, function(x) -mean(x[, "z.r"], na.rm = T))
plot(sc, meanNRI, pch = 21, bg = "forestgreen", main = "Herb: richness")
model <- lm(meanNRI ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model)

#---#

meanNRI <- sapply(ta_sc_ses, function(x) -mean(x[, "z.is"], na.ism = T))
plot(sc, meanNRI, pch = 21, bg = "tomato", main = "Tree: ind. swap")
model <- lm(meanNRI ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model)

meanNRI <- sapply(sa_sc_ses, function(x) -mean(x[, "z.is"], na.ism = T))
plot(sc, meanNRI, pch = 21, bg = "skyblue", main = "Shrub: ind. swap")
model <- lm(meanNRI ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model)

meanNRI <- sapply(ha_sc_ses, function(x) -mean(x[, "z.is"], na.ism = T))
plot(sc, meanNRI, pch = 21, bg = "forestgreen", main = "Herb: ind. swap")
model <- lm(meanNRI ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model)

#---#

meanNRI <- sapply(ta_sc_ses, function(x) -mean(x[, "z.ts"], na.ism = T))
plot(sc, meanNRI, pch = 21, bg = "tomato", main = "Tree: trial swap")
model <- lm(meanNRI ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model)

meanNRI <- sapply(sa_sc_ses, function(x) -mean(x[, "z.ts"], na.ism = T))

plot(sc, meanNRI, pch = 21, bg = "skyblue", main = "Shrub: trial swap")
model <- lm(meanNRI ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model)

meanNRI <- sapply(ha_sc_ses, function(x) -mean(x[, "z.ts"], na.ism = T))

plot(sc, meanNRI, pch = 21, bg = "forestgreen", main = "Herb: trial swap")
model <- lm(meanNRI ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model)

dev.off()



#-----------------------------------------------------------------------------#

png("figures/mean-nri-a-sc.png", 4000, 4000, pointsize = 75)

op <- par(mfrow = c(4,3), mar = c(3, 3, 2, 0.5))

meanNRI <- sapply(ta_sc_ses, function(x) -mean(x[, "z.a.tl"], na.rm = T))
plot(sc, meanNRI, pch = 21, bg = "tomato", main = "Tree: taxa labels")
model <- lm(meanNRI ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model)

meanNRI <- sapply(sa_sc_ses, function(x) -mean(x[, "z.a.tl"], na.rm = T))
plot(sc, meanNRI, pch = 21, bg = "skyblue", main = "Shrub: taxa labels")
model <- lm(meanNRI ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model)

meanNRI <- sapply(ha_sc_ses, function(x) -mean(x[, "z.a.tl"], na.rm = T))
plot(sc, meanNRI, pch = 21, bg = "forestgreen", main = "Herb: taxa labels")
model <- lm(meanNRI ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model)

#---#

meanNRI <- sapply(ta_sc_ses, function(x) -mean(x[, "z.a.r"], na.rm = T))
plot(sc, meanNRI, pch = 21, bg = "tomato", main = "Tree: richness")
model <- lm(meanNRI ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model)

meanNRI <- sapply(sa_sc_ses, function(x) -mean(x[, "z.a.r"], na.rm = T))
plot(sc, meanNRI, pch = 21, bg = "skyblue", main = "Shrub: richness")
model <- lm(meanNRI ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model)

meanNRI <- sapply(ha_sc_ses, function(x) -mean(x[, "z.a.r"], na.rm = T))
plot(sc, meanNRI, pch = 21, bg = "forestgreen", main = "Herb: richness")
model <- lm(meanNRI ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model)

#---#

meanNRI <- sapply(ta_sc_ses, function(x) -mean(x[, "z.a.is"], na.ism = T))
plot(sc, meanNRI, pch = 21, bg = "tomato", main = "Tree: ind. swap")
model <- lm(meanNRI ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model)

meanNRI <- sapply(sa_sc_ses, function(x) -mean(x[, "z.a.is"], na.ism = T))
plot(sc, meanNRI, pch = 21, bg = "skyblue", main = "Shrub: ind. swap")
model <- lm(meanNRI ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model)

meanNRI <- sapply(ha_sc_ses, function(x) -mean(x[, "z.a.is"], na.ism = T))
plot(sc, meanNRI, pch = 21, bg = "forestgreen", main = "Herb: ind. swap")
model <- lm(meanNRI ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model)

#---#

meanNRI <- sapply(ta_sc_ses, function(x) -mean(x[, "z.a.ts"], na.ism = T))
plot(sc, meanNRI, pch = 21, bg = "tomato", main = "Tree: trial swap")
model <- lm(meanNRI ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model)

meanNRI <- sapply(sa_sc_ses, function(x) -mean(x[, "z.a.ts"], na.ism = T))

plot(sc, meanNRI, pch = 21, bg = "skyblue", main = "Shrub: trial swap")
model <- lm(meanNRI ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model)

meanNRI <- sapply(ha_sc_ses, function(x) -mean(x[, "z.a.ts"], na.ism = T))

plot(sc, meanNRI, pch = 21, bg = "forestgreen", main = "Herb: trial swap")
model <- lm(meanNRI ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model)

dev.off()