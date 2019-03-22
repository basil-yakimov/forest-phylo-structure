load("clean.data/scaling.rda")
load("clean.data/scaling-ses.rda")

#_____________________________________________________________________


png("figures/regs-sc.png", 4000, 24000, pointsize = 75)

op <- par(mfrow = c(20,3))

for (ii in 1:20)
{
    hgt <- hgt_sc[[ii]]
    
    NRI <- -ta_sc_ses[[ii]][, "z"]
    plot(hgt, NRI,  pch = 21, bg = "tomato", xlab = "altitude")
    model <- lm(NRI ~ hgt)
    if (anova(model)[1, 5] < 0.05) abline(model)
    sp <- cor.test(hgt, NRI, method = "spearman", use = "complete")
    title(bquote("Tree " ~ .(ii) ~ ": " ~ rho == .({round(sp$estimate, 3)}) ~ "; " ~ p == .({round(sp$p.value, 3)}) ~ "; " ~ R^2 == .({round(cor(hgt, NRI, use = "complete")^2, 3)})))

    NRI <- -sa_sc_ses[[ii]][, "z"]
    plot(hgt, NRI,  pch = 21, bg = "tomato", xlab = "altitude")
    model <- lm(NRI ~ hgt)
    if (anova(model)[1, 5] < 0.05) abline(model)
    sp <- cor.test(hgt, NRI, method = "spearman", use = "complete")
    title(bquote("Shrub " ~ .(ii) ~ ": " ~ rho == .({round(sp$estimate, 3)}) ~ "; " ~ p == .({round(sp$p.value, 3)}) ~ "; " ~ R^2 == .({round(cor(hgt, NRI, use = "complete")^2, 3)})))
    
    NRI <- -ha_sc_ses[[ii]][, "z"]
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
  
  NRI <- -ta_sc_ses[[ii]][, "z.a"]
  plot(hgt, NRI,  pch = 21, bg = "tomato", xlab = "altitude")
  model <- lm(NRI ~ hgt)
  if (anova(model)[1, 5] < 0.05) abline(model)
  sp <- cor.test(hgt, NRI, method = "spearman", use = "complete")
  title(bquote("Tree " ~ .(ii) ~ ": " ~ rho == .({round(sp$estimate, 3)}) ~ "; " ~ p == .({round(sp$p.value, 3)}) ~ "; " ~ R^2 == .({round(cor(hgt, NRI, use = "complete")^2, 3)})))
  
  NRI <- -sa_sc_ses[[ii]][, "z"]
  plot(hgt, NRI,  pch = 21, bg = "tomato", xlab = "altitude")
  model <- lm(NRI ~ hgt)
  if (anova(model)[1, 5] < 0.05) abline(model)
  sp <- cor.test(hgt, NRI, method = "spearman", use = "complete")
  title(bquote("Shrub " ~ .(ii) ~ ": " ~ rho == .({round(sp$estimate, 3)}) ~ "; " ~ p == .({round(sp$p.value, 3)}) ~ "; " ~ R^2 == .({round(cor(hgt, NRI, use = "complete")^2, 3)})))
  
  NRI <- -ha_sc_ses[[ii]][, "z"]
  plot(hgt, NRI,  pch = 21, bg = "tomato", xlab = "altitude")
  model <- lm(NRI ~ hgt)
  if (anova(model)[1, 5] < 0.05) abline(model)
  sp <- cor.test(hgt, NRI, method = "spearman", use = "complete")
  title(bquote("Herb " ~ .(ii) ~ ": " ~ rho == .({round(sp$estimate, 3)}) ~ "; " ~ p == .({round(sp$p.value, 3)}) ~ "; " ~ R^2 == .({round(cor(hgt, NRI, use = "complete")^2, 3)})))
}

dev.off()


#_____________________________________________________________________

sc <- 1:20

# png("figures/coef-sc.png", 3000, 4000, pointsize = 75)
# 
# op <- par(mfcol = c(3,2))
# 
# tt <- c("Tree layer", "Shrub layer", "Herb layer")
# 
# for (i in 1:3){
#   scaling <- 1:length(cor_coef[i, ])
#   plot(scaling, cor_coef[i, ], pch = 21, bg = "green3", ylab = "Correlation coefficient")
#   title(tt[i])
# }
# 
# for (j in 1:3){
#   scaling <- 1:length(det_coef[j, ])
#   plot(scaling, det_coef[j, ], pch = 21, bg = "blue3", ylab = "Coefficient of determination")
#   title(tt[j])
# }
# 
# 
# dev.off()

#_____________________________________________________________________

png("figures/mean-nri-sc.png", 4000, 2500, pointsize = 75)

op <- par(mfcol = c(2,3))

meanNRI <- sapply(ta_sc_ses, function(x) -mean(x[, "z"], na.rm = T))
meanNRIa <- sapply(ta_sc_ses, function(x) -mean(x[, "z.a"], na.rm = T))

plot(sc, meanNRI, pch = 21, bg = "tomato", main = "Tree layer")
model <- lm(meanNRI ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model)

plot(sc, meanNRIa, pch = 21, bg = "tomato")
model <- lm(meanNRIa ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model)


meanNRI <- sapply(sa_sc_ses, function(x) -mean(x[, "z"], na.rm = T))
meanNRIa <- sapply(sa_sc_ses, function(x) -mean(x[, "z.a"], na.rm = T))

plot(sc, meanNRI, pch = 21, bg = "skyblue", main = "Shrub layer")
model <- lm(meanNRI ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model)

plot(sc, meanNRIa, pch = 21, bg = "skyblue")
model <- lm(meanNRIa ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model)


meanNRI <- sapply(ha_sc_ses, function(x) -mean(x[, "z"], na.rm = T))
meanNRIa <- sapply(ha_sc_ses, function(x) -mean(x[, "z.a"], na.rm = T))

plot(sc, meanNRI, pch = 21, bg = "forestgreen", main = "Herb layer")
model <- lm(meanNRI ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model)

plot(sc, meanNRIa, pch = 21, bg = "forestgreen")
model <- lm(meanNRIa ~ sc)
if (anova(model)[1, 5] < 0.05) abline(model)

dev.off()