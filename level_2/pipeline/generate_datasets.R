### R script that simulates networks under the birth-death-hybridization
### model implemented in SiPhyNetworks
### For each network, it simulates gene trees using ms
### GAB, CA, CSL (August 2022)


### Parameters ----------------------------------------------------------
#numbsim <- 150             ## number of networks to simulate
numbsim <- 5             ## number of networks to simulate
#ntips <- c(15, 30, 50)     ## number of leaves in the network
ntips <- c(15, 30)     ## number of leaves in the network
lambda <- 0.9              ## speciation rate
mu <- 0                    ## extinction rate
nu <- c(0.02, 0.04)        ## hybridization rate
hybprops <- c(1, 1, 1)     ## probabilities for each type of hybridization
#ngt <- c(100, 1000, 10000) ## number of gene trees to simulate per network
ngt <- c(5, 10, 15) ## number of gene trees to simulate per network
#gt_replics <- 30           ## number of replicates per simulating scenario
gt_replics <- 5           ## number of replicates per simulating scenario
globalseed <- 2022         ## global seed
### ---------------------------------------------------------------------

set.seed(globalseed)  
library(SiPhyNetwork) # library to simulate Networks, this use "ape" as dependence
source("functions.R") # load functions that operate on networks


# GAB: Create a directory for storing the networks and their subdirs.
# The .. is necessary as it assumes that the script generate_datasets.R
# is run from the directory pipeline.
if (!file.exists("../data")){
  dir.create(path = "../data")
}
setwd("../data")

## for testing:
#i = ntips[1]
#j = nu[1]
#k = ngt[1]

### Simulate networks with ntips leaves, nu hybridization rate and ngt gene trees
for (i in ntips) {
  for (j in nu) {
    for (k in ngt) {
      
      ## Generating the specific seed for this run:
      r_seed <- sample.int(n = 1e6, size = 1)
      cat("r_seed = ", r_seed, "\n", sep = "")
      set.seed(r_seed)

      ## Simulating networks with SiPhyNetwork
      networks <- sim.bdh.taxa.ssa(n = i,
                                   numbsim = numbsim,
                                   lambda = lambda,
                                   mu = mu ,
                                   nu = j,
                                   hybprops = hybprops,
                                   hyb.inher.fxn = make.beta.draw(1, 1),
                                   frac = 1,
                                   mrca = FALSE,
                                   complete = TRUE,
                                   stochsampling = FALSE,
                                   hyb.rate.fxn = NULL,
                                   trait.model = NULL)
      
      # GAB: Code for removing bad networks before writing
      # Get rid of null trees which go (extinct=0) and no extinct tips are sampled=1
      networks <- networks[!sapply(X = networks, FUN = is.null)]
      networks <- networks[sapply(X = networks, FUN = is.phylo)]

      # Select only networks (omit trees)
      file_networks <- vector(length = length(networks))
      for(x in 1:length(networks)){
        file_networks[x] <- as.logical(nrow(networks[[x]]$reticulation))
      }
      networks <- networks[file_networks]

      ## For every network, we simulate gene trees using PhyloCoalSimulations
      net_counter <- 1
      for (y in networks) {
        ## We do several replicates per network:
        for (ii in 1:gt_replics){
          # Folder name where all files corresponding to this network will be stored:
          filename <- paste("network", net_counter,
                          "_ntips_", i,
                          "_nu_", j,
                          "_ngt_", k, 
                          "_rep_", ii, sep="")
        
          if (!file.exists(filename)){
            dir.create(filename)
          }
          setwd(filename)
        
          # CA: Not level one warning ########
          #res <- sibCross(Tree2=y)$infor
          #ro <- any(res%in%"not level-1")
          #if(ro){
          #  sink("not_level1_warn.txt")
          #  cat(res)
          #  sink()
          #}
          ################################
        
          # generate the seed for simulate_gts.jl which uses Random.seed!(gt_seed)
          gt_seed <- sample.int(n = 1e6, size = 1)
        
          # write parenthetical format to file
          extnewick_filename <- paste(filename, ".extnewick", sep = "")
          output_newick_filename <- paste(filename, ".newick", sep = "")
          SiPhyNetwork::write.net(net = y, file = extnewick_filename)
        
          # store parenthetical format in string variable
          network_i <- write.net(net = y)

          ## writing everything in a logfile inside the folder
          logfile<-file("logfile.txt")
        
          str <- paste("ntips=",i,
            " ,lambda=",lambda,   
            " ,mu=",mu,    
            " ,nu=",j,
            " ,hybprops=",hybprops[1],",",hybprops[2],",",hybprops[3],
            " ,ngt=",k,
            ", global seed=", globalseed, 
            ", SiPhyNetwork seed=", r_seed,
            ", this is network ", net_counter, 
            ", phylocoalsims seed=", gt_seed, "\n")
        
          writeLines(str, logfile)
          #close(logfile)

          # simuate gene trees with PhyloCoalSimulations. note k = number of gene trees
          simgtcmd <- paste("julia ../../pipeline/simulate_gts.jl", extnewick_filename, output_newick_filename, k, gt_seed, sep = " ")
          system(simgtcmd)
          
          net_counter <- net_counter + 1
          setwd("../")
        }
      }
    }
  }
}


