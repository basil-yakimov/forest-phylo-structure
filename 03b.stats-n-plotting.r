load("clean.data/phylo-ses-mntd.rda")
hgt <- read.table("raw.data/hgt.txt")[[1]]

#_____________________________________________________________________


pdf("figures/tree-regs-nti.pdf", 7, 22)

op <- par(mfcol = c(7,2))

for (i in 2:ncol(t.ses.z)){
  NTI <- -t.ses.z[, i]
  plot(hgt, NTI,  pch = 21, bg = "tomato", xlab = "altitude")
  
  sp <- cor.test(hgt, NTI, method = "spearman", use = "complete")
  
  sig <- t.ses.p[, i] < 0.05 | t.ses.p[, i] > 0.95
  points(hgt[sig], NTI[sig],  pch = 21, bg = "darkred")
  
  if (sp$p.value < 0.05)
  {
    model <- lm(NTI ~ hgt)
    abline(model)
  }
  text("topleft", labels = "txt")
  legend("topleft", legend = bquote(rho == .({round(sp$estimate, 3)}) ~ "; " ~ p == .({round(sp$p.value, 3)})), bty = "n")
  
  title(main = paste0(names(t.ses.z)[i]))
}

par(op)

dev.off()


#_____________________________________________________________________


pdf("figures/shrub-regs-nti.pdf", 7, 22)

op <- par(mfcol = c(7,2))

for (i in 2:ncol(t.ses.z)){
  NTI <- -s.ses.z[, i]
  plot(hgt, NTI,  pch = 21, bg = "skyblue", xlab = "altitude")
  
  sp <- cor.test(hgt, NTI, method = "spearman", use = "complete")
  
  sig <- s.ses.p[, i] < 0.05 | s.ses.p[, i] > 0.95
  points(hgt[sig], NTI[sig],  pch = 21, bg = "darkblue")
  
  if (sp$p.value < 0.05)
  {
    model <- lm(NTI ~ hgt)
    abline(model)
  }
  text("topleft", labels = "txt")
  legend("topleft", legend = bquote(rho == .({round(sp$estimate, 3)}) ~ "; " ~ p == .({round(sp$p.value, 3)})), bty = "n")
  
  title(main = paste0(names(s.ses.z)[i]))
}

par(op)

dev.off()

#_____________________________________________________________________


pdf("figures/herb-regs-nti.pdf", 7, 22)

op <- par(mfcol = c(7,2))

for (i in 2:ncol(t.ses.z)){
  NTI <- -h.ses.z[, i]
  plot(hgt, NTI,  pch = 21, bg = "limegreen", xlab = "altitude")
  
  sp <- cor.test(hgt, NTI, method = "spearman", use = "complete")
  
  sig <- h.ses.p[, i] < 0.05 | h.ses.p[, i] > 0.95
  points(hgt[sig], NTI[sig],  pch = 21, bg = "darkgreen")
  
  if (sp$p.value < 0.05)
  {
    model <- lm(NTI ~ hgt)
    abline(model)
  }
  text("topleft", labels = "txt")
  legend("topleft", legend = bquote(rho == .({round(sp$estimate, 3)}) ~ "; " ~ p == .({round(sp$p.value, 3)})), bty = "n")
  
  title(main = paste0(names(h.ses.z)[i]))
}

par(op)

dev.off()
