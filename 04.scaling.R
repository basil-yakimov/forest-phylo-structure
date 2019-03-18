load("clean.data/tree-abund.rda")
load("clean.data/shrub-abund.rda")
load("clean.data/herb-abund.rda")

hgt <- read.table("raw.data/hgt.txt")[[1]]


a <- vector(mode = "list", length = 3)
a[[1]] <- ta
a[[2]] <- sa
a[[3]] <- ha


#________________________________________

#considering the boundaries between transects

#1110-1180,1190-1255, 1260-1338, 1344-1408, 1414-1496, 1504-1552, 1558-1616, 1622-1682, and 1687-1734

#pairs <- matrix(c(11, 12, 20, 21, 33, 34, 44, 45, 58, 59, 66, 67, 76, 77, 87, 88), ncol = 2, byrow = T)

#________________________________________

out_tot <- vector(mode = "list", length = 3)

for (i in 1:3){
  dt <- a[[i]]
  
  sc <- 1:20
  
  out2 <- vector(mode = "list", length = length(sc))
  hgt_sc <- vector(mode = "list", length = 7) 
  
  for (jj in 1:length(sc))
  {
    dd <- matrix(0, nrow = (96-sc[jj]+1), ncol = ncol(dt))
    h <- NULL
    del <- NULL
    for (ii in 1:(96-sc[jj]+1))
    {
      dd[ii, ] <- colSums(dt[(ii):(ii+sc[jj]-1), ])
      del <- c(del, ii %in% seq(1, 96, by = sc[jj]))
    }
    
    for (ii in seq(1, 96-sc[jj]+1, by = sc[jj])){
      h <- c(h, mean(hgt[(ii):(ii+sc[jj]-1)]))
    }
    
    dd <- dd[del, ]
    
    hgt_sc[[jj]] <- h
    
    colnames(dd) <- colnames(dt)
    out2[[jj]] <- dd
  }
  
  out_tot[[i]] <- out2
}

save(out_tot, file = "clean.data/scaling-abund.rda")
save(hgt_sc, file = "clean.data/scaling-hgt.rda")
