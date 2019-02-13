load("clean.data/scaling-2-phylo-ses.rda")
hgt <- read.table("raw.data/hgt.txt")[[1]]
count <- c(2:11, 13:20, 22:33, 35:44, 46:58, 60:66, 68:76, 78:87, 89:96)
hgt2 <- c()
for (ii in count){
  hgt2 <- c(hgt2, mean(hgt[ii-1], hgt[ii]))
}
hgt <- hgt2

#_____________________________________________________________________


pdf("figures/tree-regs-2.pdf", 7, 22)

op <- par(mfcol = c(7,2))

for (i in 2:ncol(t.ses.z2)){
  NRI <- -t.ses.z2[, i]
  plot(hgt, NRI,  pch = 21, bg = "tomato", xlab = "altitude")
  
  sp <- cor.test(hgt, NRI, method = "spearman", use = "complete")
  
  sig <- t.ses.p2[, i] < 0.05 | t.ses.p2[, i] > 0.95
  points(hgt[sig], NRI[sig],  pch = 21, bg = "darkred")
  
  if (sp$p.value < 0.05)
  {
    model <- lm(NRI ~ hgt)
    abline(model)
  }
  text("topleft", labels = "txt")
  legend("topleft", legend = bquote(rho == .({round(sp$estimate, 3)}) ~ "; " ~ p == .({round(sp$p.value, 3)})), bty = "n")
  
  title(main = paste0(names(t.ses.z2)[i]))
}

par(op)

dev.off()


#_____________________________________________________________________


pdf("figures/shrub-regs-2.pdf", 7, 22)

op <- par(mfcol = c(7,2))

for (i in 2:ncol(t.ses.z2)){
  NRI <- -s.ses.z2[, i]
  plot(hgt, NRI,  pch = 21, bg = "skyblue", xlab = "altitude")
  
  sp <- cor.test(hgt, NRI, method = "spearman", use = "complete")
  
  sig <- s.ses.p2[, i] < 0.05 | s.ses.p2[, i] > 0.95
  points(hgt[sig], NRI[sig],  pch = 21, bg = "darkblue")
  
  if (sp$p.value < 0.05)
  {
    model <- lm(NRI ~ hgt)
    abline(model)
  }
  text("topleft", labels = "txt")
  legend("topleft", legend = bquote(rho == .({round(sp$estimate, 3)}) ~ "; " ~ p == .({round(sp$p.value, 3)})), bty = "n")
  
  title(main = paste0(names(s.ses.z2)[i]))
}

par(op)

dev.off()

#_____________________________________________________________________


pdf("figures/herb-regs-2.pdf", 7, 22)

op <- par(mfcol = c(7,2))

for (i in 2:ncol(t.ses.z2)){
  NRI <- -h.ses.z2[, i]
  plot(hgt, NRI,  pch = 21, bg = "limegreen", xlab = "altitude")
  
  sp <- cor.test(hgt, NRI, method = "spearman", use = "complete")
  
  sig <- h.ses.p2[, i] < 0.05 | h.ses.p2[, i] > 0.95
  points(hgt[sig], NRI[sig],  pch = 21, bg = "darkgreen")
  
  if (sp$p.value < 0.05)
  {
    model <- lm(NRI ~ hgt)
    abline(model)
  }
  text("topleft", labels = "txt")
  legend("topleft", legend = bquote(rho == .({round(sp$estimate, 3)}) ~ "; " ~ p == .({round(sp$p.value, 3)})), bty = "n")
  
  title(main = paste0(names(h.ses.z2)[i]))
}

par(op)

dev.off()
