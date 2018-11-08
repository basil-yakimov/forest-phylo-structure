load("workspaces/shrub-ses.RData")

#_____________________________________________________________________

s.coef.pe <- c()
s.coef.sp <- c()
s.coef.ke <- c()
s.coef.p.pe <- c()
s.coef.p.sp <- c()
s.coef.p.ke <- c()

pdf("output/shrubs-coef.pdf", 22, 7)

old.par <- par(no.readonly = TRUE)
par(mfrow = c(2,7))

for (i in 2:ncol(s.ses)){
  a <- s.ses[, i]
  plot(s.ses$hgt, a,  pch = 21, bg = "steelblue")
  model <- lm(a ~ s.ses$hgt)
  abline(model)
  
  pe <- cor.test(s.ses$hgt, a)
  sp <- cor.test(s.ses$hgt, a, method = "spearman")
  ke <- cor.test(s.ses$hgt, a, method = "kendall")
  
  title(main = paste0(names(s.ses)[i], " - ",
                      round(pe$estimate, dig = 3), " - ", round(pe$p.value, 3),
                      round(sp$estimate, dig = 3), " - ", round(sp$p.value, 3),
                      round(ke$estimate, dig = 3), " - ", round(ke$p.value, 3)))
  
  s.coef.pe <- c(s.coef.pe, pe$estimate)
  s.coef.sp <- c(s.coef.sp, sp$estimate)
  s.coef.ke <- c(s.coef.ke, ke$estimate)
  s.coef.p.pe <- c(s.coef.p.pe, pe$p.value)
  s.coef.p.sp <- c(s.coef.p.sp, sp$p.value)
  s.coef.p.ke <- c(s.coef.p.ke, ke$p.value)
}

par(old.par)

dev.off()

s.coef <- data.frame(null.model = names(s.ses)[2:ncol(s.ses)], pearson.cor.coef = s.coef.pe,
                     pearson.p.value = s.coef.p.pe,
                     spearman.cor.coef = s.coef.sp,
                     spearman.p.value = s.coef.p.sp,
                     kendall.cor.coef = s.coef.ke,
                     kendall.p.value = s.coef.p.ke)

save(s.coef, file = "workspaces/shrub-coef.rda")
