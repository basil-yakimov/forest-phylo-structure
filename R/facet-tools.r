t.decomp <- function(dat, q)
{
    N <- dim(dat)[1]
    
    x <- dat[,colSums(dat) > 0]
    p <- x/rowSums(x)
    
    z <- colSums(p)/N
    #z <- z/sum(z)
    
    if (q == 0) 
    {
        Dg <- sum(z > 0)
        Da <- sum(p > 0)/N
        Cqn <-  (N - Dg/Da)/(N-1)
        Uqn <- (Da/Dg - 1/N) / (1 - 1/N)
    } else
    
    
    if (q == 1)
    {
        Dg <-  - sum( z * log(z), na.rm = T)
        Da <- - sum( p * log(p), na.rm = T) / N
        
        #w <- rowSums(x)
        #w <- w/sum(w)
        #( Da - Dg -  sum( w * log(w), na.rm = T) ) / log(N)
        
        Cqn <- Uqn <- 1 - (Dg-Da)/log(N)
        
        Dg <- exp(Dg)
        Da <- exp(Da)
    } else
    
    {
        Dg <- sum( z^q ) ^ (1/(1-q))
        Da <- ( sum( p^q )/N ) ^ (1/(1-q))
        Cqn <- (N^(1-q) - (Dg/Da)^(1-q)) / (N^(1-q) - 1)
        Uqn <- ( (Da/Dg)^(1-q) - N^(q-1) ) / ( 1 - N^(q-1) )
    }
    
    
    res <- c(a = Da, b = Dg/Da, g = Dg, c = Cqn, u = Uqn, q = q)
    
    return(res)
}


p.decomp <- function(d, tree, q)
{
    require(geiger)
    
    x <- d[,colSums(d) > 0]
    x <- x/rowSums(x)
    
    z <- colSums(x)
    N <- sum(z)
    
    if (q == 1) q <- 0.9999
    if (q == 0) q <- 10^(-10)
    
    tmp.tree <- treedata(tree, t(x), warnings = F)$phy
    T <- max(branching.times(tmp.tree))
    branches <- matrix(NA, nrow = nrow(tmp.tree$edge), ncol = 5)
    branches[,1:2] <- tmp.tree$edge
    branches[,3] <- tmp.tree$edge.length
    for (ii in 1:nrow(branches))
    {
        leaves.node = tips(tmp.tree, branches[ii, 2])
        branches[ii, 4] <- (sum(z[leaves.node]/N, na.rm = T)/T)^q
        
        if (length(leaves.node) > 1) branches[ii, 5] <- sum( rowSums(x[,leaves.node]/N/T)^q, na.rm = T) else
            branches[ii, 5] <- sum( (x[,leaves.node]/N/T)^q, na.rm = T)
    }
    
    
    
    PDg <- sum(branches[,3]*branches[,4])^(1/(1-q))
    PDa <- sum(branches[,3]*branches[,5])^(1/(1-q))/N
    
    Cqn <- (N^(1-q) - (PDg/PDa)^(1-q)) / (N^(1-q) - 1)
    Uqn <- ( (PDa/PDg)^(1-q) - N^(q-1) ) / ( 1 - N^(q-1) )


    res <- c(a = PDa, b = PDg/PDa, g = PDg, c = Cqn, u = Uqn, q = q)

    return(res)
}


t.dist <- function(d, q, index = "c")
{
    N <- dim(d)[1]
    res <- matrix(rep(0, N*N), nrow = N)
    rownames(res) <- colnames(res) <- rownames(d)
    for (ii in 1:N)
    {
        for (jj in 1:N)
        {
            res[ii,jj] <- t.decomp(d[c(ii,jj),], q)[index]
        }
    }
    
    return(as.dist(1-res))
}

p.dist <- function(d, tree, q, index = "c")
{
    N <- dim(d)[1]
    res <- matrix(rep(0, N*N), nrow = N)
    rownames(res) <- colnames(res) <- rownames(d)
    for (ii in 1:N)
    {
        for (jj in 1:N)
        {
            res[ii,jj] <- p.decomp(d[c(ii,jj),], tree, q)[index]
        }
    }
    
    return(as.dist(1-res))
}
