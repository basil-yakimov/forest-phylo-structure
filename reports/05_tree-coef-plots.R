load("workspaces/tree-ses.RData")

#tses.tl.m <- lm(taxa.labels ~ hgt, t.ses)
#plot(t.ses$hgt, t.ses$taxa.labels, pch = 21, bg = "steelblue")
#abline(tses.tl.m)
#title(main = paste0("taxa.labels - ", round(pe$estimate, dig = 3), " - ", round(pe$p.value, 3),
#                    round(sp$estimate, dig = 3), " - ", round(sp$p.value, 3),
#                    round(ke$estimate, dig = 3), " - ", round(ke$p.value, 3)))

#_____________________________________________________________________

t.coef.pe <- c()
t.coef.sp <- c()
t.coef.ke <- c()
t.coef.p.pe <- c()
t.coef.p.sp <- c()
t.coef.p.ke <- c()

pdf("output/trees-coef.pdf", 22, 7)

old.par <- par(no.readonly = TRUE)
par(mfrow = c(2,7))

for (i in 2:ncol(t.ses)){
  a <- t.ses[, i]
  plot(t.ses$hgt, a,  pch = 21, bg = "steelblue")
  model <- lm(a ~ t.ses$hgt)
  abline(model)
  
  pe <- cor.test(t.ses$hgt[!is.na(a)], a[!is.na(a)])
  sp <- cor.test(t.ses$hgt[!is.na(a)], a[!is.na(a)], method = "spearman")
  ke <- cor.test(t.ses$hgt[!is.na(a)], a[!is.na(a)], method = "kendall")
  
  title(main = paste0(names(t.ses)[i], " - ",
                      round(pe$estimate, dig = 3), " - ", round(pe$p.value, 3),
                      round(sp$estimate, dig = 3), " - ", round(sp$p.value, 3),
                      round(ke$estimate, dig = 3), " - ", round(ke$p.value, 3)))
  
  t.coef.pe <- c(t.coef.pe, pe$estimate)
  t.coef.sp <- c(t.coef.sp, sp$estimate)
  t.coef.ke <- c(t.coef.ke, ke$estimate)
  t.coef.p.pe <- c(t.coef.p.pe, pe$p.value)
  t.coef.p.sp <- c(t.coef.p.sp, sp$p.value)
  t.coef.p.ke <- c(t.coef.p.ke, ke$p.value)
}

par(old.par)

dev.off()

t.coef <- data.frame(null.model = names(t.ses)[2:ncol(t.ses)], pearson.cor.coef = t.coef.pe,
                     pearson.p.value = t.coef.p.pe,
                     spearman.cor.coef = t.coef.sp,
                     spearman.p.value = t.coef.p.sp,
                     kendall.cor.coef = t.coef.ke,
                     kendall.p.value = t.coef.p.ke)

save(t.coef, file = "workspaces/tree-coef.rda")
