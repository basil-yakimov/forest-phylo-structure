library(phytools)

my.add <- function (tree, species, genus = NULL, where = c("root", "random"), tol) 
{
  if (!inherits(tree, "phylo")) 
    stop("tree should be an object of class \"phylo\".")
  if (!is.ultrametric(tree, tol = tol)) 
    warning("this code has only been tested with ultrametric tree\n  your tree may be returned without edge lengths")
  where <- where[1]
  if (is.null(genus)) {
    x <- strsplit(species, "")[[1]]
    i <- 1
    while (x[i] != "_" && x[i] != " ") i <- i + 1
    genus <- paste(x[2:i - 1], collapse = "")
  }
  ii <- grep(paste(genus, "_", sep = ""), tree$tip.label)
  if (length(ii) > 1) {
    if (!is.monophyletic(tree, tree$tip.label[ii])) 
      warning(paste(genus, "may not be monophyletic\n  attaching to the most inclusive group containing members of this genus"))
    nn <- findMRCA(tree, tree$tip.label[ii])
    if (where == "root") 
      tree <- my.bind.tip(tree, gsub(" ", "_", species), tol = tol,
                       where = nn)
    else if (where == "random") {
      tt <- splitTree(tree, list(node = nn, bp = tree$edge.length[which(tree$edge[, 
                                                                                  2] == nn)]))
      tt[[2]] <- add.random(tt[[2]], tips = gsub(" ", 
                                                 "_", species))
      tree <- paste.tree(tt[[1]], tt[[2]])
    }
    else stop("option 'where' not recognized")
  }
  else if (length(ii) == 1) {
    nn <- ii
    if (where == "root") 
      tree <- my.bind.tip(tree, gsub(" ", "_", species), tol = tol,
                       where = nn, position = 0.5 * tree$edge.length[which(tree$edge[, 
                                                                                     2] == nn)])
    else if (where == "random") 
      tree <- bind.tip(tree, gsub(" ", "_", species), tol = tol,
                       where = nn, position = runif(n = 1) * tree$edge.length[which(tree$edge[, 
                                                                                              2] == nn)])
    else stop("option 'where' not recognized")
  }
  else warning("could not match your species to a genus\n  check spelling, including case")
  tree
}


my.bind.tip <- function (tree, tip.label, edge.length = NULL, where = NULL, tol,
                         position = 0, interactive = FALSE, ...) 
{
  if (!inherits(tree, "phylo")) 
    stop("tree should be an object of class \"phylo\".")
  use.edge.length <- if (is.null(tree$edge.length)) 
    FALSE
  else TRUE
  if (use.edge.length == FALSE) 
    tree <- compute.brlen(tree)
  if (interactive == TRUE) {
    plotTree(tree, ...)
    cat(paste("Click where you would like to bind the tip \"", 
              tip.label, "\"\n", sep = ""))
    flush.console()
    obj <- get.treepos(message = FALSE)
    where <- obj$where
    position <- obj$pos
  }
  else if (is.null(where)) 
    where <- length(tree$tip.label) + 1
  if (where <= length(tree$tip.label) && position == 0) {
    pp <- 1e-12
    if (tree$edge.length[which(tree$edge[, 2] == where)] <= 
        1e-12) {
      tree$edge.length[which(tree$edge[, 2] == where)] <- 2e-12
      ff <- TRUE
    }
    else ff <- FALSE
  }
  else pp <- position
  if (is.null(edge.length) && is.ultrametric(tree, tol = tol)) {
    H <- nodeHeights(tree)
    if (where == (length(tree$tip.label) + 1)) 
      edge.length <- max(H)
    else edge.length <- max(H) - H[tree$edge[, 2] == where, 
                                   2] + position
  }
  tip <- list(edge = matrix(c(2, 1), 1, 2), tip.label = tip.label, 
              edge.length = edge.length, Nnode = 1)
  class(tip) <- "phylo"
  obj <- bind.tree(tree, tip, where = where, position = pp)
  if (where <= length(tree$tip.label) && position == 0) {
    nn <- obj$edge[which(obj$edge[, 2] == which(obj$tip.label == 
                                                  tip$tip.label)), 1]
    obj$edge.length[which(obj$edge[, 2] == nn)] <- obj$edge.length[which(obj$edge[, 
                                                                                  2] == nn)] + 1e-12
    obj$edge.length[which(obj$edge[, 2] == which(obj$tip.label == 
                                                   tip$tip.label))] <- 0
    obj$edge.length[which(obj$edge[, 2] == which(obj$tip.label == 
                                                   tree$tip.label[where]))] <- 0
  }
  root.time <- if (!is.null(obj$root.time)) 
    obj$root.time
  else NULL
  obj <- untangle(obj, "read.tree")
  if (!is.null(root.time)) 
    obj$root.time <- root.time
  if (interactive) 
    plotTree(obj, ...)
  if (!use.edge.length) 
    obj$edge.length <- NULL
  obj
}
