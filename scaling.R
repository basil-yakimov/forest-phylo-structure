load("clean.data/tree-abund.rda")
load("clean.data/shrub-abund.rda")
load("clean.data/herb-abund.rda")

hgt <- read.table("raw.data/hgt.txt")[[1]]


a <- vector(mode = "list", length = 3)
a[[1]] <- ta
a[[2]] <- sa
a[[3]] <- ha


pairs <- matrix(c(11, 12, 56, 57), ncol = 2, byrow = T)

sc <- 2:15

out <- vector(mode = "list", length = length(sc))

for (jj in 1:length(sc))
{
  dd <- matrix(0, nrow = 96-sc[jj]+1, ncol = 149)
  del <- NULL
  for (ii in 1:(96-sc[jj]+1))
  {
    dd[ii, ] <- colSums(ha[(ii):(ii+sc[jj]-1), ])
    for (uu in 1:nrow(pairs))
    {
      if ( pairs[uu, 1]  %in% ((ii):(ii+sc[jj]-1)) & pairs[uu, 2] %in% (ii):(ii+sc[jj]-1) )
      {
        del <- c(del, ii)
      }
    }
  }
  
  del <- unique(del)
  
  out[[jj]] <- dd[-del, ]
}



#________________________________________

#considering the boundaries between transects

#1110-1180,1190-1255, 1260-1338, 1344-1408, 1414-1496, 1504-1552, 1558-1616, 1622-1682, and 1687-1734

#hgt <- read.table("raw.data/hgt.txt")[[1]]
#h <- row.names(ta)
#names(h) <- hgt
pairs <- matrix(c(11, 12, 20, 21, 33, 34, 44, 45, 58, 59, 66, 67, 76, 77, 87, 88), ncol = 2, byrow = T)

h <- row.names(ta)
names(h) <- hgt #looking on the numbers of sites
count <- c(2:11, 13:20, 22:33, 35:44, 46:58, 60:66, 68:76, 78:87, 89:96)

ta2 <- ta
for (ii in count){
  a <- ta[ii-1, ] + ta[ii, ]
  for (i in 1:length(a)){   #filling the row
    ta2[ii, i] <- a[i]
  }
}
ta2 <- ta2[count, ] #removing residues


sc1 <- 2:8

out1_tot <- vector(mode = "list", length = 3)

for (i in 1:3){
  dt <- a[[i]]
  out <- vector(mode = "list", length = length(sc1))
  hgt_sc <- vector(mode = "list", length = 7) 
  
  for (jj in 1:length(sc1))
  {
    dd <- matrix(0, nrow = 96-sc1[jj]+1, ncol = ncol(dt)) 
    del <- NULL 
    h <- NULL
    
    for (ii in 1:(96-sc1[jj]+1))
    {
      dd[ii, ] <- colSums(dt[(ii):(ii+sc1[jj]-1), ]) 
      h <- c(h, mean(hgt[(ii):(ii+sc1[jj]-1)]))
      
      for (uu in 1:nrow(pairs))
      {
        if ( pairs[uu, 1]  %in% ((ii):(ii+sc1[jj]-1)) & pairs[uu, 2] %in% (ii):(ii+sc1[jj]-1) )
        {
          del <- c(del, ii)
        }
      }
    }
    
    del <- unique(del)
    h <- h[-del]
    hgt_sc[[jj]] <- h
    
    colnames(dd) <- colnames(dt)
    out[[jj]] <- dd[-del, ] 
    
  }
  
  out1_tot[[i]] <- out
}



#________________________________________

#ignoring the boundaries between transects


out2_tot <- vector(mode = "list", length = 3)

for (i in 1:3){
  dt <- a[[i]]
  
  sc <- 2:nrow(dt)
  
  out2 <- vector(mode = "list", length = length(sc))
  hgt_sc2 <- vector(mode = "list", length = 7) 
  
  for (jj in 1:length(sc))
  {
    dd <- matrix(0, nrow = 96-sc[jj]+1, ncol = ncol(dt))
    h <- NULL
    for (ii in 1:(96-sc[jj]+1))
    {
      dd[ii, ] <- colSums(dt[(ii):(ii+sc[jj]-1), ])
      h <- c(h, mean(hgt[(ii):(ii+sc[jj]-1)]))
    }
    
    hgt_sc2[[jj]] <- h
    
    colnames(dd) <- colnames(dt)
    out2[[jj]] <- dd
  }
  
  out2_tot[[i]] <- out2
}




#________________________________________


save(out1_tot, out2_tot, file = "clean.data/scaling-abund.rda")
save(hgt_sc, hgt_sc2, file = "clean.data/scaling-hgt.rda")
