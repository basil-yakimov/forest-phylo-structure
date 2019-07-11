plot.ses <- function(ses, hgt = hgt, col = "grey", lab = "", xlab = "hgt") {
  plot(hgt, ses, pch = 19, col = col, xlab = xlab, ylab = lab)
  fit1 <- lm(ses ~ hgt)
  fit2 <- lm(ses ~ hgt + I(hgt^2))
  sm1 <- summary(fit1)
  sm2 <- summary(fit2)
  p1 <- pf(sm1$fstatistic[1], df1 = sm1$fstatistic[2], df2 = sm1$fstatistic[3], lower.tail = F)
  p2 <- pf(sm2$fstatistic[1], df1 = sm2$fstatistic[2], df2 = sm2$fstatistic[3], lower.tail = F)
  usr <- par("usr")
  if (p1 < 0.05)
  {
    if (p2 >= 0.05)
    {
      abline(fit1)
      text(usr[1], usr[3] + (usr[4] - usr[3]) / 10, round(sm1$r.squared, digits = 3), pos = 4, adj = c(1, 1))
    }
    else
    {
      n <- length(ses)
      delta <- AIC(fit2) + (2*4*5/(n-4-1)) - AIC(fit1) - (2*3*4/(n-3-1))
      if (delta < 0)
      {
        x <- seq(min(hgt), max(hgt), len = 1000)
        abc <- coef(fit2)
        y <- abc[1] + abc[2] * x + abc[3] * x^2
        lines(x, y)
        text(usr[1], usr[3] + (usr[4] - usr[3]) / 10, round(sm2$r.squared, digits = 3), pos = 4, adj = c(1, 1))
      }
      else
      {
        abline(fit1)
        text(usr[1], usr[3] + (usr[4] - usr[3]) / 10, round(sm1$r.squared, digits = 3), pos = 4, adj = c(1, 1))
      }
    }
  }
  else if (p2 < 0.05)
  {
    x <- seq(min(hgt), max(hgt), len = 1000)
    abc <- coef(fit2)
    y <- abc[1] + abc[2] * x + abc[3] * x^2
    lines(x, y)
    text(usr[1], usr[3] + (usr[4] - usr[3]) / 10, round(sm2$r.squared, digits = 3), pos = 4, adj = c(1, 1))
  }
  
}