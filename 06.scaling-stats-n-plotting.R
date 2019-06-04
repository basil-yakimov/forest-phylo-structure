load("clean.data/scaling.rda")
load("clean.data/scaling-ses.rda")

#_____________________________________________________________________



plot.ses <- function(ses, hgt = hgt, col, lab) {
  plot(hgt, ses, pch = 19, col = col, ylab = lab)
  fit1 <- lm(ses ~ hgt)
  fit2 <- lm(ses ~ hgt + I(hgt^2))
  sm1 <- summary(fit1)
  sm2 <- summary(fit2)
  p1 <- pf(sm1$fstatistic[1], df1 = sm1$fstatistic[2], df2 = sm1$fstatistic[3], lower.tail = F)
  p2 <- pf(sm2$fstatistic[1], df1 = sm2$fstatistic[2], df2 = sm2$fstatistic[3], lower.tail = F)
  usr <- par("usr")
  if (p1 < 0.05)
  {
    if (p2 >= 0.05)
    {
      abline(fit1)
      text(usr[1], usr[3] + (usr[4] - usr[3]) / 10, round(sm1$r.squared, digits = 3), pos = 4, adj = c(1, 1))
    }
    else
    {
      n <- length(ses)
      delta <- AIC(fit2) + (2*4*5/(n-4-1)) - AIC(fit1) - (2*3*4/(n-3-1))
      if (delta < 0)
      {
        x <- seq(min(hgt), max(hgt), len = 1000)
        abc <- coef(fit2)
        y <- abc[1] + abc[2] * x + abc[3] * x^2
        lines(x, y)
        text(usr[1], usr[3] + (usr[4] - usr[3]) / 10, round(sm2$r.squared, digits = 3), pos = 4, adj = c(1, 1))
      }
      else
      {
        abline(fit1)
        text(usr[1], usr[3] + (usr[4] - usr[3]) / 10, round(sm1$r.squared, digits = 3), pos = 4, adj = c(1, 1))
      }
    }
  }
  else if (p2 < 0.05)
  {
    x <- seq(min(hgt), max(hgt), len = 1000)
    abc <- coef(fit2)
    y <- abc[1] + abc[2] * x + abc[3] * x^2
    lines(x, y)
    text(usr[1], usr[3] + (usr[4] - usr[3]) / 10, round(sm2$r.squared, digits = 3), pos = 4, adj = c(1, 1))
  }
  
}

png("figures/regs-sc.png", 4000, 24000, pointsize = 75)

op <- par(mfrow = c(20,3))

for (ii in 1:20)
{
  hgt <- hgt_sc[[ii]]
  NRI <- -ta_sc_ses[[ii]][, "z.is"]
  
  plot.ses(ses = NRI, hgt = hgt, col = "tomato", lab = "NRI")
  
  title(bquote("Tree " ~ .(ii)))
  
  NRI <- -sa_sc_ses[[ii]][, "z.is"]
  
  plot.ses(ses = NRI, hgt = hgt, col = "skyblue", lab = "NRI")
  
  title(bquote("Shrub " ~ .(ii)))
  
  NRI <- -ha_sc_ses[[ii]][, "z.is"]
  
  plot.ses(ses = NRI, hgt = hgt, col = "forestgreen", lab = "NRI")
  
  title(bquote("Herb " ~ .(ii)))
}

dev.off()


png("figures/regs-sc-a.png", 4000, 24000, pointsize = 75)

op <- par(mfrow = c(20,3))

for (ii in 1:20)
{
  hgt <- hgt_sc[[ii]]
  NRI <- -ta_sc_ses[[ii]][, "z.a.is"]
  
  plot.ses(ses = NRI, hgt = hgt, col = "tomato", lab = "NRI")
  
  title(bquote("Tree " ~ .(ii)))
  
  NRI <- -sa_sc_ses[[ii]][, "z.a.is"]
  
  plot.ses(ses = NRI, hgt = hgt, col = "skyblue", lab = "NRI")
  
  title(bquote("Shrub " ~ .(ii)))
  
  NRI <- -ha_sc_ses[[ii]][, "z.a.is"]
  
  plot.ses(ses = NRI, hgt = hgt, col = "forestgreen", lab = "NRI")
  
  title(bquote("Herb " ~ .(ii)))
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



#__________________________________________________________________________________

