### R script that simulates networks under the birth-death-hybridization
### model implemented in SiPhyNetworks
### For each network, it simulates gene trees using ms
### GAB, CA, CSL (August 2022)

## toy=true uses small number of networks and gene trees for debugging
toy=TRUE

### Parameters ----------------------------------------------------------
lambda <- 0.9              ## speciation rate
mu <- 0                    ## extinction rate
nu <- c(0.02, 0.04)        ## hybridization rate
hybprops <- c(1, 1, 1)     ## probabilities for each type of hybridization
globalseed <- 2022         ## global seed

if(toy){
  numbsim <- 5             ## number of networks to simulate
  ntips <- c(15, 30)       ## number of leaves in the network
  ngt <- c(5, 10)      ## number of gene trees to simulate per network
  gt_replics <- 5          ## number of replicates per simulating scenario
}else{
  numbsim <- 150             ## number of networks to simulate
  ntips <- c(15, 30, 50)     ## number of leaves in the network
  ngt <- c(100, 1000, 10000) ## number of gene trees to simulate per network
  gt_replics <- 30           ## number of replicates per simulating scenario
}
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
                                   mrca = TRUE,
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
          filename <- paste("net", net_counter,
                          "_ntips_", i,
                          "_nu_", j,
                          "_ngt_", k, 
                          "_rep_", ii, sep="")
        
          if (!file.exists(filename)){
            dir.create(filename)
          }
          setwd(filename)
                
          # generate the seed for simulate_gts.jl which uses Random.seed!(gt_seed)
          gt_seed <- sample.int(n = 1e6, size = 1)
          writeLines(as.character(gt_seed), "gt_seed")
          
          # write number of gene trees to file for sim_gts_calc_cf_startingtree.jl
          writeLines(as.character(k), "num_gt")
        
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
          close(logfile)          
          net_counter <- net_counter + 1
          setwd("../")
        }
      }
    }
  }
}
setwd("../pipeline")
system("julia sim_gts_calc_cf_startingtree.jl ../data")



# Creating directory to save the results compresed to send to chtc
dir.create("../compess")
dir.create("../compess/files_15_tips")
dir.create("../compess/files_20_tips")
dir.create("../compess/files_30_tips")
files <- dir("../data")
nam <- "job"
file_nam <- "file"

pos_15 <- regexpr(pattern="ntips_15",text=dir("../data"))
files_15 <- files[pos_15!=-1]
for(i in 1:length(files_15)){
  #i=1
  dir_compress <-paste("../data/", files_15[i],sep="")
  dir_save <- paste("../compess/files_15_tips/",nam,i,sep="")
  dir.create(dir_save)
  setwd(dir_compress)
  run <- paste("tar -czvf"," ../",dir_save ,"/",file_nam,".tar.gz",
               " `ls *.extnewick *.csv *.tre gt_seed logfile.txt`",sep="")
  system(run)
  
  setwd("../../pipeline")
}

pos_20 <- regexpr(pattern="ntips_20",text=dir("../data"))
files_20 <- files[pos_20!=-1]
for(i in 1:length(files_20)){
  #i=1
  dire <-paste("../data/", files_20[i],sep="")
  lugar <- paste("../compess/files_20_tips/",nam,i,sep="")
  dir.create(lugar)
  setwd(dire)
  run <- paste("tar -czvf"," ../",lugar,"/file.tar.gz"," `ls`",sep="")
  system(run)
  
  setwd("../../pipeline")
}

pos_30 <- regexpr(pattern="ntips_30",text=dir("../data"))
files_30 <- files[pos_30!=-1]
for(i in 1:length(files_30)){
  #i=1
  dire <-paste("../data/", files_30[i],sep="")
  lugar <- paste("../compess/files_30_tips/",nam,i,sep="")
  dir.create(lugar)
  setwd(dire)
  run <- paste("tar -czvf"," ../",lugar,"/file.tar.gz"," `ls`",sep="")
  system(run)
  
  setwd("../../pipeline")
}
