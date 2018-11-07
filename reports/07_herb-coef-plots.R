load("workspaces/herb-ses.rda")

#_____________________________________________________________________

h.coef.pe <- c()
h.coef.sp <- c()
h.coef.ke <- c()
h.coef.p.pe <- c()
h.coef.p.sp <- c()
h.coef.p.ke <- c()

pdf("output/herb-coef.pdf", 22, 7)

old.par <- par(no.readonly = TRUE)
par(mfrow = c(2,7))

for (i in 2:ncol(h.ses)){
  a <- h.ses[, i]
  plot(h.ses$hgt, a,  pch = 21, bg = "steelblue")
  model <- lm(a ~ h.ses$hgt)
  abline(model)
  
  pe <- cor.test(h.ses$hgt, a)
  sp <- cor.test(h.ses$hgt, a, method = "spearman")
  ke <- cor.test(h.ses$hgt, a, method = "kendall")
  
  title(main = paste0(names(h.ses)[i], " - ",
                      round(pe$estimate, dig = 3), " - ", round(pe$p.value, 3),
                      round(sp$estimate, dig = 3), " - ", round(sp$p.value, 3),
                      round(ke$estimate, dig = 3), " - ", round(ke$p.value, 3)))
  
  h.coef.pe <- c(h.coef.pe, pe$estimate)
  h.coef.sp <- c(h.coef.sp, sp$estimate)
  h.coef.ke <- c(h.coef.ke, ke$estimate)
  h.coef.p.pe <- c(h.coef.p.pe, pe$p.value)
  h.coef.p.sp <- c(h.coef.p.sp, sp$p.value)
  h.coef.p.ke <- c(h.coef.p.ke, ke$p.value)
}

par(old.par)

dev.off()

h.coef <- data.frame(null.model = names(h.ses)[2:ncol(h.ses)], pearson.cor.coef = h.coef.pe,
                     pearson.p.value = h.coef.p.pe,
                     spearman.cor.coef = h.coef.sp,
                     spearman.p.value = h.coef.p.sp,
                     kendall.cor.coef = h.coef.ke,
                     kendall.p.value = h.coef.p.ke)

save(h.coef, file = "workspaces/herb-coef.rda")