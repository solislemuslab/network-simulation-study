####GAB Parameters were specified on notion in https://www.notion.so/Project-on-SNaQ-limitations-80e419d6d58a42308e4bf0f5cceb9752#e5df0d1a4c97433ca14f2e6388be966a
numbsim <- 150
ntips <- c(15, 30, 50)
lambda <- 0.9
mu <- 0
nu <- c(0.02, 0.04)
hybprops <- c(1, 1, 1) # we are starting with all-equal
ngt <- c(100, 1000, 10000)
set.seed(2022) # our initial seed, all the others below will depend on sample.int()

#GAB  original construction of the is.phylo function in geiger. Unnecessary to load the whole thing for just a function
is.phylo <- function(x) {
    "phylo" %in% class(x)
}

library(SiPhyNetwork) #library to simulate Networks, this use "ape" as dependence

#GAB create a directory for storing the networks and their subdirs. the .. is necessary as it assumes that the script generate_datasets.R us run from the directory pipeline, which should only contain code and not data
dir.create(path = "../data")


#3.1 Simulate networks
#GAB set each seed as below, They are all pseudo-random, but dependent on the first set.seed(2022) above

#GAB nested for for generating the params one at a time for those which vary, i.e., ntips, nu, ngt

for (i in ntips) {
    for (j in nu) {
        for (k in ngt) {
            setwd("../data")
            r_seed <- sample.int(n = 1e6, size = 1)
            cat("r_seed = ", r_seed, "\n", sep = "")
            set.seed(r_seed)

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
            #GAB Code for removing bad networks before writing
            # get rid of null trees which go extinct=0 and no extinct tips are sampled=1
            networks <- networks[!sapply(X = networks, FUN = is.null)]
            networks <- networks[sapply(X = networks, FUN = is.phylo)]
                        
            #Code to select only networks (omit trees)
            file_networks <- vector(length = length(networks))
            for(x in 1:length(networks)){          
                file_networks[x] <- as.logical(nrow(networks[[x]]$reticulation))
            }
            networks <- networks[file_networks]
            
            net_counter <- 1
            for (y in networks) {
                filename <- paste("network", net_counter,
                                  "_ntips_", i,
                                  "_nu_", j,
                                  "_ngt_", k, sep="")

                dir.create(filename)
                setwd(filename)
                
                gt_seed <- sample.int(n = 1e6, size = 1)
                cat("gt_seed = ", gt_seed, "\n", sep = "")
                extnewick_filename <- paste(filename, ".extnewick", sep = "")
                
                # write to string file
                SiPhyNetwork::write.net(net = y, file = extnewick_filename)
                
                hybridlambda_filename <- paste(filename, ".hybridlambda", sep = "")
                
                # convert from extnewick to hybridlambda
                system(paste("julia extnewick2hybridlambda.jl ", extnewick_filename, " ", hybridlambda_filename, sep = ""))

                # run hybrid-Lambda on filename_extnewick and capture the output
                system(paste("hybrid-Lambda -spcu ", "'",
                             net_hyla_format,"'", " -num ", 
                             k," -seed ", gt_seed, " -o ",
                             hybridlambda_filename," > hybridlambda_output 2>&1",sep=""))

                notultram_bool <- sum(grepl(x = readLines("hybridlambda_output"),
                                        pattern = "WARNING! NOT ULTRAMETRIC!!!"))
                notultram_bool
                
                if(notultram_bool){#Is ultrametric
				  #GAB the network file in extnewick will be called appending the parameter values as well as a counter for numbering each network from 1 no length(networks)
                    setwd("../")
                    system(paste("rm -rf ", filename, sep = ""))
                } else {
                    file.remove("hybridlambda_output")
                    setwd("../")
                }
                net_counter <- net_counter + 1
            }
        }
    }
}

