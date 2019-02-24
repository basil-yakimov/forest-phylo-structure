load("clean.data/scaling-phylo-ses.rda")
load("clean.data/scaling-hgt.rda")

#_____________________________________________________________________


png("figures/tree-regs-sc.png", 4000, 8000, pointsize = 75)

op <- par(mfcol = c(7,3))

for (jj in 1:3){
  for (ii in 1:7){
    NRI <- -c[[jj]][[ii]][, 1]
    
    hgt <- hgt_sc[[ii]]
    
    plot(hgt, NRI,  pch = 21, bg = "tomato", xlab = "altitude")
    
    sp <- cor.test(hgt, NRI, method = "spearman", use = "complete")
    
    if (sp$p.value < 0.05)
    {
      model <- lm(NRI ~ hgt)
      abline(model)
    }
    text("topleft", labels = "txt")
    legend("topleft", legend = bquote(rho == .({round(sp$estimate, 3)}) ~ "; " ~ p == .({round(sp$p.value, 3)})), bty = "n")
  }
}

dev.off()
