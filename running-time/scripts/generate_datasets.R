library(SiPhyNetwork)
library(ape)

### Revised simulation strategy as per 2022 (from notion)
#GAB 1. R script with simulation of networks and gene trees
#GAB     1. Input: Parameters (lambda, mu, hybprobs, N, ntips, ngt, nu, seed)
#GAB     2. Code:
#GAB         2.1. Simulate trees with ape
#GAB         2.2. For every tree, simulate gene trees with HybridLambda
#GAB     3. Output: One folder (name with parameter values) with N subfolders (one per simulated tree). Each subfolder has the simulated tree (.hybridlambda) and another text file with the simulated gene trees (.genetrees)
#GAB 2. Julia script with snaq runs
#GAB     1. Input: One subfolder from before (with one simulated network and its simulated gene trees)
#GAB     2. Code:
#GAB         1. Generate CF table
#GAB         2. Run snaq for h=0 nruns=10 (parallel)
#GAB         3. Run snaq for h=1,2,3,...until either 5 or the max number of hybridizations in the simulated network
#GAB     3. Output: One folder (name with parameter values) with subfolders (one per h value in snaq) with the output files from snaq

#####################################################################################
#1. Input Parameters for simulating nerworks and gene trees
#####################################################################################
numbsim <- 150
ntips <- c(15, 30, 50)
lambda <- 0.9
mu <- 0
nu <- c(0.02, 0.04)
ngt <- c(100, 1000, 10000)
set.seed(2022) # our initial seed, all the others below will depend on sample.int()

setwd("../data/simulations")
#GAB set each seed as below, They are all pseudo-random, but dependent on the first set.seed(2022) above
set.seed(sample.int(n = 1e6, size = 1))

#GAB nested for for generating the params one at a time for those which vary, i.e., ntips, nu, ngt

for (i in ntips) {
    for (j in nu) {
        for (k in ngt) {
            # CODE FOR SIOMULATINF COALESCENT TREES, rtree(type = "coal") or rcoal() or something from ape. change nested structure accordingly
            #GAB Code for removing bad trees before writing
            # get rid of null trees which go extinct=0 and no extinct tips are sampled=1
            trees <- bdNS1[!sapply(X = trees, FUN = is.null)]
            trees <- bdNS1[sapply(X = trees, FUN = is.phylo)]
            #GAB code for picking only the trees of interest, Level-2? Ultrametric? We are now interested in taking only the good nets onwards
            # CODE HERE

            #GAB after selecting only the trees that we need, run the julia script extnewick2hybridlambda.jl for format conversion over each element tree in trees
            net_counter <- 1
            for (tree in trees) {
                #GAB the tree file in extnewick will be called appending the parameter values as well as a counter for numbering each tree from 1 no length(trees)
                tree_filename <- paste("tree", tree_counter, "_ntips", i, "_nu", j, "_ngt", k)
                # create the subdirectory for a single tree using tree_filename
                dir.create(tree_filename)
                tree_counter <- tree_counter + 1
                # write the tree to a file in extnewick format, inside the directory tree_filename
                ape::write.[tree CHECK](net = network, file = paste(tree_filename, "/", tree_filename, ".newick", sep =""))
                # run the script for format conversion in julia using the extnewick filename and an output filename as input
                # INPUT would be something like net1_ntips15_nu0.2_ngt100/net1_ntips15_nu0.2_ngt100.extnewick 
                # OUTPUT would be something like net1_ntips15_nu0.2_ngt100/net1_ntips15_nu0.2_ngt100.hybridlambda
                system(command = paste("julia extnewick2hybridlambda.jl ", INPUT, " ", OUTPUT, sep = ""))
            }
        }
    }
}

#GAB Up to here, we have a series of directories netX_ntipsI_nuJ_ngtK, each with two files, one .extnewick, one .hybridlambda. Now we will iterate over each directory, and run hybrid-Lambda over each .hybridlambda file
for i in (dir(pattern = "tree")) {
    # visit the dir i
    setwd(i)
    # recreate the filename this time only with the .hybridlambda one
    tree_fielname <- dir(pattern = ".hybridlambda")
    # run hybrid-Lambda on filename
    command <- paste("hybrid-Lambda -spcu ", tree_filename,
                     " -num ", ngt,
                     " -seed ", sample.int(n = 1e6, size = 1),
                     " -o ", tree_filename, ".genetrees",
                     sep="")
    system(command)
    # get back to the parent directory
    setwd("..")
}
