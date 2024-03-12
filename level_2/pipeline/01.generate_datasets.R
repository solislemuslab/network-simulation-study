### R script that simulates networks under the birth-death-hybridization
### model implemented in SiPhyNetworks using maximum granularity

### Original Parameters ----------------------------------------------------------
#lambda <- 0.9              ## speciation rate
#mu <- 0                    ## extinction rate
#nu <- c(0.02, 0.04)        ## hybridization rate
#hybprops <- c(1, 1, 1)     ## probabilities for each type of hybridization THIS NEEDS TO BE PROVIDED AS A STRING SEPARATED BY COMMAS
#globalseed <- 2022         ## global seed
#numbsim <- 150             ## number of networks to simulate
#ntips <- c(15, 30, 50)     ## number of leaves in the network
#ngt <- c(100, 1000, 10000) ## number of gene trees to simulate per network
#gt_replics <- 30           ## number of replicates per simulating scenario

### ---------------------------------------------------------------------
# ordered arguments for command line:
#ntips
#lambda
#mu
#nu
#hybprops
#seed
#outpath

args <- commandArgs(TRUE)

ntips <- as.numeric(args[1])
lambda <- as.numeric(args[2])
mu <- as.numeric(args[3])
nu <- as.numeric(args[4])
hybprops <- as.numeric(unlist(strsplit(args[5], split=",")))
seed <- as.numeric(args[6])
outpath = args[7]

library(SiPhyNetwork) # library to simulate Networks, this use "ape" as dependence
# starting execution in scripts
source("functions.R") # load functions that operate on networks

# set output dir
setwd(outpath)

# set the seed
set.seed(seed)

# set the while conditional for starting the simulation
continue <- TRUE
maxiter <- 1000
counter <- 1

while (continue & (counter <= maxiter)) {
    network <- sim.bdh.taxa.ssa(n = ntips,
                                numbsim = 1,
                                lambda = lambda,
                                mu = mu ,
                                nu = nu,
                                hybprops = hybprops,
                                hyb.inher.fxn = make.beta.draw(1, 1),
                                frac = 1,
                                twolineages = TRUE,
                                complete = TRUE,
                                stochsampling = FALSE,
                                hyb.rate.fxn = NULL,
                                trait.model = NULL)[[1]]
    # check whether the network is phylo or restart if not
    isphylo <- is.phylo(network)
    if (!isphylo) {
        counter <- counter + 1
        next
    }
    # now check whether it is network or restart if not
    isnetwork <- as.logical(nrow(network$reticulation))
    if (!isnetwork) {
        counter <- counter + 1
        next
    }
    # if it is, write to file and break the while
    if (isnetwork) {
        SiPhyNetwork::write.net(net = network, file = "network.extnewick")
        cat("Attempted ", counter, " times until successfully picking a network\n", sep = "")
        break
    }
}
if (counter > maxiter) {
    stop("Simulation was unsuccessful with ", maxiter, " attepmts:\n  Try rising maxiter\n", sep="")
}
