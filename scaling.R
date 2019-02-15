load("clean.data/tree-abund.rda")
load("clean.data/shrub-abund.rda")
load("clean.data/herb-abund.rda")

hgt <- read.table("raw.data/hgt.txt")[[1]]


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

#1110-1180,1190-1255, 1260-1338, 1344-1408, 1414-1496, 1504-1552, 1558-1616, 1622-1682, and 1687-1734
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

#________________________________________

sa2 <- sa
for (ii in count){
  a <- sa[ii-1, ] + sa[ii, ]
  for (i in 1:length(a)){
    sa2[ii, i] <- a[i]
  }
}
sa2 <- sa2[count, ]

#________________________________________

ha2 <- ha
for (ii in count){
  a <- ha[ii-1, ] + ha[ii, ]
  for (i in 1:length(a)){
    ha2[ii, i] <- a[i]
  }
}
ha2 <- ha2[count, ]

#________________________________________

ha <- ha2
sa <- sa2
ta <- ta2

save(ta, sa, ha, file = "clean.data/scaling-2-abund.rda")
