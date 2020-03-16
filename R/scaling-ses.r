library(picante)

rand.mpd <- function(cdm, dist, ab = F)
{
  perm <- sample(nrow(dist))
  names(cdm) <- names(cdm)[perm]
  mpd(cdm, dist, ab)
}


scaling.grain <- function(cdm, dis, ab = F, iter = 100, silent = F)
{
  K <- nrow(cdm)
  
  rmax <- floor((K-1)/2)
  n <- sum(K - 2*(1:rmax))
  
  res <- data.frame(matrix(NA, nrow = n, ncol = 3))
  colnames(res) <- c("r", "id", "ses")
  
  counter <- 1
  
  for (r in 1:rmax)
  {
    nloc <- K - 2*r
    for (jj in 1:nloc)
    {
      if (!silent) print(paste0("r = ", r, ", num = ", jj))
      
      res$r[counter]  <- r
      res$id[counter] <- r + jj
      
      target <- as.data.frame(t(colSums(cdm[jj : (jj + 2*r), ])))
      
      null.model <- replicate(iter, rand.mpd(target, dis, ab))
      res$ses[counter] <- (mpd(target, dis, ab) - mean(null.model))/sd(null.model)
      counter <- counter + 1
    }
  }
  return(res)
}


scaling.pool <- function(cdm, dis, ab = F, iter = 100, silent = F)
{
  K <- nrow(cdm)
  
  rmax <- floor((K-1)/2)
  n <- sum(K - 2*(1:rmax))
  
  res <- data.frame(matrix(NA, nrow = n, ncol = 3))
  colnames(res) <- c("r", "id", "ses")
  
  counter <- 1
  
  for (r in 1:rmax)
  {
    nloc <- K - 2*r
    for (jj in 1:nloc)
    {
      if (!silent) print(paste0("r = ", r, ", num = ", jj))
      
      res$r[counter]  <- r
      res$id[counter] <- r + jj
      
      target <- cdm[r + jj, ]
      pool <- colSums(cdm[jj : (jj + 2*r), ]) > 0

      if (sum(pool) > 1)
      {
        target <- target[pool]
        pool.dist <- dis[pool, pool]
        
        null.model <- replicate(iter, rand.mpd(target, pool.dist, ab))
        res$ses[counter] <- (mpd(target, pool.dist, ab) - mean(null.model))/sd(null.model)
      }

      counter <- counter + 1
    }
  }
  return(res)
}


scaling.coherent <- function(cdm, dis, ab = F, iter = 100, silent = F)
{
  K <- nrow(cdm)
  
  rmax <- floor((K-1)/6)
  n <- sum(K - 2*(1:rmax) - 2*(2*(1:rmax) + 1))
  
  res <- data.frame(matrix(NA, nrow = n, ncol = 3))
  colnames(res) <- c("r", "id", "ses")
  
  counter <- 1
  
  for (r in 1:rmax)
  {
    nloc <- K - 2*r - 2*(2*r + 1)
    for (jj in 1:nloc)
    {
      if (!silent) print(paste0("r = ", r, ", num = ", jj))
      
      res$r[counter]  <- r
      res$id[counter] <- r + jj + 2*r + 1
      
      target <- as.data.frame(t(colSums(cdm[(2*r + 1 + jj) : (2*r + 1 + jj + 2*r), ])))
      pool <- colSums(cdm[jj : (jj + 6*r + 2), ]) > 0
      
      if (sum(pool) > 1)
      {
        target <- target[pool]
        pool.dist <- dis[pool, pool]
        
        null.model <- replicate(iter, rand.mpd(target, pool.dist, ab))
        res$ses[counter] <- (mpd(target, pool.dist, ab) - mean(null.model))/sd(null.model)
      }
      counter <- counter + 1
    }
  }
  return(res)
}


plot.scaling <- function(frame, col = "wheat")
{
  require(nlme)
  require(effects)
  
  assign(".frame", frame, env = .GlobalEnv)
  
  fit <- lme(fixed = ses ~ r, random = ~ 1 | id, data = .frame, method = "REML", na.action = "na.omit")
  ef <- effect("r", fit)
  plot(ses ~ r, .frame, col = col, pch = 19, xlab = "Масштаб, размер блока", ylab = "SES MPD")
  for (ii in 1:max(.frame$id)) lines(ses ~ r, .frame, subset = id == ii, col = "grey")
  if (anova(fit)[2, 4] < 0.05)
  {
    lines(ef$x$r, ef$fit, type = "l", lwd = 3, ylim = range(ef$lower, ef$upper), pch = 16)
    lines(ef$x$r, ef$lower, lwd = 2, col = "grey40")
    lines(ef$x$r, ef$upper, lwd = 2, col = "grey40")
  }
  legend("topright", legend = "", bty = "n", title = paste0("coef: ", round(coef(fit)[1, 2], 3), "; p-value: ", round(anova(fit)[2, 4], 3)))
}
