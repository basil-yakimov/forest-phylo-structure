load("clean.data/scaling-phylo-ses.rda")
load("clean.data/scaling-hgt.rda")

#_____________________________________________________________________


png("figures/tree-regs-sc.png", 4000, 24000, pointsize = 75)

op <- par(mfcol = c(20,3))

cor_coef <- matrix(0, nrow = 3, ncol = 20)
det_coef <- matrix(0, nrow = 3, ncol = 20)

for (jj in 1:3){
  for (ii in 1:20){
    NRI <- -s[[jj]][[ii]][, 1]
    
    hgt <- hgt_sc[[ii]]
    
    plot(hgt, NRI,  pch = 21, bg = "tomato", xlab = "altitude")
    
    sp <- cor.test(hgt, NRI, method = "spearman", use = "complete")
    model <- lm(NRI ~ hgt)
    
    if (s[[jj]][[ii]][, 2] < 0.05)
    {
      abline(model)
    }
    
    text("topleft", labels = "txt")
    legend("topleft", legend = bquote(rho == .({round(sp$estimate, 3)}) ~ "; " ~ p == .({round(sp$p.value, 3)})), bty = "n")
    
    cor_coef[jj, ii] <- cor(hgt, NRI, method = "spearman")
    det_coef[jj, ii] <- cor(hgt, NRI)^2
    legend("topright", legend = bquote("R^2" == .({round(cor(hgt, NRI)^2, 3)})), bty = "n")
    
    if (jj == 1) {title(c("Tree layer ", ii))}
    if (jj == 2) {title(c("Shrub layer ", ii))}
    if (jj == 3) {title(c("Herb layer ", ii))}
  }
}

dev.off()

#_____________________________________________________________________

png("figures/coef-sc.png", 3000, 4000, pointsize = 75)

op <- par(mfcol = c(3,2))

tt <- c("Tree layer", "Shrub layer", "Herb layer")

for (i in 1:3){
  scaling <- 1:length(cor_coef[i, ])
  plot(scaling, cor_coef[i, ], pch = 21, bg = "green3", ylab = "Correlation coefficient")
  title(tt[i])
}

for (j in 1:3){
  scaling <- 1:length(det_coef[j, ])
  plot(scaling, det_coef[j, ], pch = 21, bg = "blue3", ylab = "Coefficient of determination")
  title(tt[j])
}


dev.off()

#_____________________________________________________________________

png("figures/mean-nri-sc.png", 4000, 3500, pointsize = 75)

tt <- c("Tree layer", "Shrub layer", "Herb layer")

op <- par(mfcol = c(3,3))

for (jj in 1:3){
  meanNRI <- NULL
  meanNRI_low <- NULL
  meanNRI_high <- NULL
  
  for (ii in 1:20){
   meanNRI <- c(meanNRI, mean(-s[[jj]][[ii]][, 1]))
   meanNRI_low <- c(meanNRI_low, mean(-s[[jj]][[ii]][1:(length(s[[jj]][[ii]][, 1])%/%2), 1]))
   meanNRI_high <- c(meanNRI_high, mean(-s[[jj]][[ii]][(1+length(s[[jj]][[ii]][, 1])%/%2):(length(s[[jj]][[ii]][, 1])), 1]))
  }
  
  scaling <- 1:length(meanNRI)
  plot(scaling, meanNRI, pch = 21, bg = rgb(0.7, (jj/4), 0.2, 1))
  title(tt[jj])
  
  model <- lm(meanNRI ~ scaling)
  abline(model)
  
  plot(scaling, meanNRI_low, pch = 21, bg = rgb(0.7, (jj/4), 0.2, 1))
  plot(scaling, meanNRI_high, pch = 21, bg = rgb(0.7, (jj/4), 0.2, 1))
  
}

dev.off()

#_____________________________________________________________________

rm(list = ls())

load("clean.data/scaling-phylo-ses-abund.rda")

png("figures/mean-nri-a-sc.png", 4000, 3500, pointsize = 75)

tt <- c("Tree layer", "Shrub layer", "Herb layer")

op <- par(mfcol = c(3,3))

for (jj in 1:3){
  meanNRI <- NULL
  meanNRI_low <- NULL
  meanNRI_high <- NULL
  
  for (ii in 1:20){
    meanNRI <- c(meanNRI, mean(-s[[jj]][[ii]][, 1]))
    meanNRI_low <- c(meanNRI_low, mean(-s[[jj]][[ii]][1:(length(s[[jj]][[ii]][, 1])%/%2), 1]))
    meanNRI_high <- c(meanNRI_high, mean(-s[[jj]][[ii]][(1+length(s[[jj]][[ii]][, 1])%/%2):(length(s[[jj]][[ii]][, 1])), 1]))
  }
  
  scaling <- 1:length(meanNRI)
  plot(scaling, meanNRI, pch = 21, bg = rgb(0.7, (jj/4), 0.2, 1))
  title(tt[jj])
  
  model <- lm(meanNRI ~ scaling)
  abline(model)
  
  plot(scaling, meanNRI_low, pch = 21, bg = rgb(0.7, (jj/4), 0.2, 1))
  plot(scaling, meanNRI_high, pch = 21, bg = rgb(0.7, (jj/4), 0.2, 1))
  
}

dev.off()
