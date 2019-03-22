load("clean.data/tree-abund.rda")
load("clean.data/shrub-abund.rda")
load("clean.data/herb-abund.rda")

hgt <- read.table("raw.data/hgt.txt")[[1]]

#________________________________________

#considering the boundaries between transects

#1110-1180,1190-1255, 1260-1338, 1344-1408, 1414-1496, 1504-1552, 1558-1616, 1622-1682, and 1687-1734

#pairs <- matrix(c(11, 12, 20, 21, 33, 34, 44, 45, 58, 59, 66, 67, 76, 77, 87, 88), ncol = 2, byrow = T)

#________________________________________

sc <- 1:20

hgt_sc <- ta_sc <- sa_sc <- ha_sc <- vector(mode = "list", length = length(sc)) 

for (jj in 1:length(sc))
{
  num <- 96 %/% sc[jj]
  
  dt <- matrix(0, nrow = num, ncol = ncol(ta))
  ds <- matrix(0, nrow = num, ncol = ncol(sa))
  dh <- matrix(0, nrow = num, ncol = ncol(ha))
  dhgt <- rep(NA, num)
  
  colnames(dt) <- colnames(ta)
  colnames(ds) <- colnames(sa)
  colnames(dh) <- colnames(ha)
  
  for (ii in 1:num)
  {
    dt[ii, ] <- colSums(ta[((ii-1)*sc[jj]+1):(ii*sc[jj]), ])
    ds[ii, ] <- colSums(sa[((ii-1)*sc[jj]+1):(ii*sc[jj]), ])
    dh[ii, ] <- colSums(ha[((ii-1)*sc[jj]+1):(ii*sc[jj]), ])
    dhgt[ii] <- mean(hgt[((ii-1)*sc[jj]+1):(ii*sc[jj])])
  }
  
  ta_sc[[jj]] <- dt
  sa_sc[[jj]] <- ds
  ha_sc[[jj]] <- dh
  hgt_sc[[jj]] <- dhgt
}

save(hgt_sc, ta_sc, sa_sc, ha_sc, file = "clean.data/scaling.rda")