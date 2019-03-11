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
