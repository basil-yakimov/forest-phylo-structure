load("clean.data/phylo-ses-mntd.rda")
hgt <- read.table("raw.data/hgt.txt")[[1]]

#_____________________________________________________________________


pdf("figures/tree-regs-2.pdf", 7, 22)

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
    legend("topright", legend = round(coef(model)[2]*100, 4), bty = "n")
  }
  text("topleft", labels = "txt")
  legend("topleft", legend = bquote(rho == .({round(sp$estimate, 3)}) ~ "; " ~ p == .({round(sp$p.value, 3)})), bty = "n")
  
  title(main = paste0(names(t.ses.z)[i]))
}

par(op)

dev.off()


#_____________________________________________________________________


pdf("figures/shrub-regs-2.pdf", 7, 22)

op <- par(mfcol = c(7,2))

for (i in 2:ncol(t.ses.z)){
  NTI <- -s.ses.z[, i]
  plot(hgt, NTI,  pch = 21, bg = "skyblue", xlab = "altitude")
  
  sp <- cor.test(hgt, NTI, method = "spearman", use = "complete")
  
  sig <- s.ses.p[, i] < 0.05 | s.ses.p[, i] > 0.95
  points(hgt[sig], NTI[sig],  pch = 21, bg = "darkblue")
  
  model <- lm(NTI ~ hgt)
  if (sp$p.value < 0.05)
  {
    abline(model)
  }
  legend("topright", legend = round(coef(model)[2]*100, 4), bty = "n")
  
  text("topleft", labels = "txt")
  legend("topleft", legend = bquote(rho == .({round(sp$estimate, 3)}) ~ "; " ~ p == .({round(sp$p.value, 3)})), bty = "n")
  
  title(main = paste0(names(s.ses.z)[i]))
}

par(op)

dev.off()

#_____________________________________________________________________


pdf("figures/herb-regs-2.pdf", 7, 22)

op <- par(mfcol = c(7,2))

for (i in 2:ncol(t.ses.z)){
  NTI <- -h.ses.z[, i]
  plot(hgt, NTI,  pch = 21, bg = "limegreen", xlab = "altitude")
  
  sp <- cor.test(hgt, NTI, method = "spearman", use = "complete")
  
  sig <- h.ses.p[, i] < 0.05 | h.ses.p[, i] > 0.95
  points(hgt[sig], NTI[sig],  pch = 21, bg = "darkgreen")
  
  model <- lm(NTI ~ hgt)
  if (sp$p.value < 0.05)
  {
    abline(model)
  }
  legend("topright", legend = round(coef(model)[2]*100, 4), bty = "n")
  
  text("topleft", labels = "txt")
  legend("topleft", legend = bquote(rho == .({round(sp$estimate, 3)}) ~ "; " ~ p == .({round(sp$p.value, 3)})), bty = "n")
  
  title(main = paste0(names(h.ses.z)[i]))
}

par(op)

dev.off()


#_____________________________________________________________________
#

png("figures/tree-nti.png", 600, 400)

par(mar = c(4.5, 4.5, 1.5, 1))

NRI <- -t.ses.z[, 2]
plot(hgt, NRI, xlab = "Высота над уровнем моря, м", cex = 2, cex.axis = 2, cex.lab = 2, pch = 21, bg = "tomato")
abline(lm(NRI ~ hgt), lwd = 2)
points(hgt[t.ses.p[, 2] < 0.05 | t.ses.p[, 2] > 0.95],
       NRI[t.ses.p[, 2] < 0.05 | t.ses.p[, 2] > 0.95],
       cex = 2, pch = 21, bg = "darkred")

dev.off()


png("figures/shrub-nti.png", 600, 400)

par(mar = c(4.5, 4.5, 1.5, 1))

NRI <- -s.ses.z[, 2]
plot(hgt, NRI, xlab = "Высота над уровнем моря, м", cex = 2, cex.axis = 2, cex.lab = 2, pch = 21, bg = "skyblue")
points(hgt[s.ses.p[, 2] < 0.05 | s.ses.p[, 2] > 0.95],
       NRI[s.ses.p[, 2] < 0.05 | s.ses.p[, 2] > 0.95],
       cex = 2, pch = 21, bg = "darkblue")

dev.off()


png("figures/herb-nti.png", 600, 400)

par(mar = c(4.5, 4.5, 1.5, 1))

NRI <- -h.ses.z[, 2]
plot(hgt, NRI, xlab = "Высота над уровнем моря, м", cex = 2, cex.axis = 2, cex.lab = 2, pch = 21, bg = "limegreen")
points(hgt[h.ses.p[, 2] < 0.05 | h.ses.p[, 2] > 0.95],
       NRI[h.ses.p[, 2] < 0.05 | h.ses.p[, 2] > 0.95],
       cex = 2, pch = 21, bg = "darkgreen")

dev.off()


#_____________________________________________________________________
#

output <- matrix(0, nrow = 14, ncol = 3)
colnames(output) <- c("spearman", "p-value", "r^2")

for (i in 2:ncol(t.ses.z)){
  NRI <- -t.ses.z[, i]
  
  output[i-1, 1] <- round(cor.test(hgt, NRI, method = "spearman")$estimate, 3)
  output[i-1, 2] <- round(cor.test(hgt, NRI, method = "spearman", use = "complete")$p.value, 3)
  output[i-1, 3] <- round(cor.test(hgt, NRI)$estimate^2, 3)
  
}

rownames(output) <- colnames(t.ses.z)[2:15]
write.csv2(output, file = "figures/tree-nti.csv")



output <- matrix(0, nrow = 14, ncol = 3)
colnames(output) <- c("spearman", "p-value", "r^2")

for (i in 2:ncol(s.ses.z)){
  NRI <- -s.ses.z[, i]
  
  output[i-1, 1] <- round(cor.test(hgt, NRI, method = "spearman")$estimate, 3)
  output[i-1, 2] <- round(cor.test(hgt, NRI, method = "spearman", use = "complete")$p.value, 3)
  output[i-1, 3] <- round(cor.test(hgt, NRI)$estimate^2, 3)
  
}

rownames(output) <- colnames(t.ses.z)[2:15]
write.csv2(output, file = "figures/shrub-nti.csv")



output <- matrix(0, nrow = 14, ncol = 3)
colnames(output) <- c("spearman", "p-value", "r^2")

for (i in 2:ncol(h.ses.z)){
  NRI <- -h.ses.z[, i]
  
  output[i-1, 1] <- round(cor.test(hgt, NRI, method = "spearman")$estimate, 3)
  output[i-1, 2] <- round(cor.test(hgt, NRI, method = "spearman", use = "complete")$p.value, 3)
  output[i-1, 3] <- round(cor.test(hgt, NRI)$estimate^2, 3)
  
}

rownames(output) <- colnames(t.ses.z)[2:15]
write.csv2(output, file = "figures/herb-nti.csv")
