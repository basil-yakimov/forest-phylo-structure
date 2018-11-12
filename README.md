# forest-phylo-structure

A set of scripts and data for the analysis of phylogenetic structure of three vegetation layers in Donglishan mountains.
We use phylogenetic tree from [Zanne et al., 2014](http://www.nature.com/nature/journal/v506/n7486/full/nature12872.html).

* First, we clean the data and construct two phylogenetic trees for woody and herbaceous species. We clean latin names for species present in Zanne phylogeny and add polytomies to genera in cases where species present in our data are absent in Zanne phylogeny.
* Second, we perform phylogenetic community structure analysis with package picante tools.
* Third, we perform statistical analisys of the relation of net relatedness index (NRI) to altitude.

The main result is that phylogenetic clustering increases with altitude. This effect is consistent between layers.