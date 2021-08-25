library(ape)
library(phytools)

# this code unfortunately assumes that simulation objects come from
# the binary_tree_simulations.R script, which saves a copy of these
# objects for the next line to load them. uncomment if needed
load("TreeSim_NetSim.Rdata")
#source("binary_tree_simulations.R")

### summarize binary trees for TreeSim in each experiment

# create an object for results of iterative testing whether trees are binary
bins_TreeSim <- list()
length(bins_TreeSim) <- length(simTreeSim)
index <- 1

# change the class for each experiment list to multiPhylo
simTreeSim <- lapply(simTreeSim, FUN = function(x) `class<-`(x, "multiPhylo"))

for (i in simTreeSim) {
    # remove elements that are not of class phylo
    i <- i[sapply(X = i, FUN = function(x) inherits(x, "phylo"))]
    # remove elements that are null
    i <- i[!sapply(X = i, FUN = is.null)]
    # append the freq table of binary to non-binary trees to output object
    bins_TreeSim[index] <- list(table(is.binary(i)))
    # update index
    index <- index + 1
}

# update names and show results for TreeSim
# please note that numbers other than 100 represent cases where
# some simulated elements were either non-phylo or null
names(bins_TreeSim) <- names(simTreeSim)
print(bins_TreeSim)
#$TS_TRUE_4
#
#TRUE 
# 100 
#
#$TS_TRUE_6
#
#TRUE 
# 100 
#
#$TS_TRUE_8
#
#TRUE 
# 100 
#
#$TS_FALSE_4
#
#TRUE 
#  98 
#
#$TS_FALSE_6
#
#TRUE 
# 100 
#
#$TS_FALSE_8
#
#TRUE 
# 100

### summarize binary trees for NetSim in each experiment

# create an object for results of iterative testing whether trees are binary
bins_NetSim <- list()
length(bins_NetSim) <- length(simNetSim)
index <- 1

simNetSim <- lapply(simNetSim, FUN = function(x) `class<-`(x, "multiPhylo"))

for (i in simNetSim) {
    # remove elements that are not of class phylo
    i <- i[sapply(X = i, FUN = function(x) inherits(x, "phylo"))]
    # remove elements that are null
    i <- i[!sapply(X = i, FUN = is.null)]
    # append the freq table of binary to non-binary trees to output object
    bins_NetSim[index] <- list(table(is.binary(i)))
    # update index
    index <- index + 1
}

# update names and show results for NetSim
# please note that numbers other than 100 represent cases where
# some simulated elements were either non-phylo or null
names(bins_NetSim) <- names(simNetSim)
print(bins_NetSim)
#$NS_TRUE_4
#
#TRUE 
# 100 
#
#$NS_TRUE_6
#
#TRUE 
# 100 
#
#$NS_TRUE_8
#
#TRUE 
# 100 
#
#$NS_FALSE_4
#
#FALSE 
#  100 
#
#$NS_FALSE_6
#
#FALSE 
#  100 
#
#$NS_FALSE_8
#
#FALSE 
#  100
