
####GAB Parameters were specified on notion in https://www.notion.so/Project-on-SNaQ-limitations-80e419d6d58a42308e4bf0f5cceb9752#e5df0d1a4c97433ca14f2e6388be966a
numbsim <- 150
#numbsim <- 50
ntips <- c(15, 30, 50)
#ntips <- c(15)
lambda <- 0.9
mu <- 0
nu <- c(0.02, 0.04)
#nu <- c(0.02)
hybprops <- c(1, 1, 1) # we are starting with all-equal
ngt <- c(100, 1000, 10000)
#ngt <- c(50)
gt_replics <- 30 # we want 30 replicates of gene tree samples
set.seed(2022) # our initial seed, all the others below will depend on sample.int()

#GAB  original construction of the is.phylo function in geiger. Unnecessary to load the whole thing for just a function
is.phylo <- function(x) {
  "phylo" %in% class(x)
}

library(SiPhyNetwork) #library to simulate Networks, this use "ape" as dependence
source("functions.R") # load functions that operate on networks


#setwd(patscripts)
#GAB create a directory for storing the networks and their subdirs. the .. is necessary as it assumes that the script generate_datasets.R us run from the directory pipeline, which should only contain code and not data
dir.create(path = "../data")


#3.1 Simulate networks
#GAB set each seed as below, They are all pseudo-random, but dependent on the first set.seed(2022) above

#GAB nested for for generating the params one at a time for those which vary, i.e., ntips, nu, ngt

for (i in ntips) {
  for (j in nu) {
    for (k in ngt) {
      #i=ntips[1];j=nu[1];k=ngt[1]
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
        
        #net_counter <- 2;y=networks[[net_counter]]
        filename <- paste("network", net_counter,
                          "_ntips_", i,
                          "_nu_", j,
                          "_ngt_", k, sep="")
        
        dir.create(filename)
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
        
        # CA We have to use it as seed to generare gene tree with "ms"
        gt_seed <- sample.int(n = 1e6, size = 1)
        
        extnewick_filename <- paste(filename, ".extnewick", sep = "")
        
        # write to string file
        SiPhyNetwork::write.net(net = y, file = extnewick_filename)
        
        #CA
        network_i <- write.net(net = y)
        
        #CA converting networks to ms format
        system(paste('ms-converter --newick ', "'", network_i , "'>ms_convert_res.txt", sep = ""))
        
        #CA reading ms format file and change the gene tree number (by default is 1)
        #CA ms 15 1 -T -I 15 1 1 1......
        #CA for do that split out form  position when -T appear to the end of the string
        #split out from the beginning the string until the number 1 that make reference to the gene tree number 
        #and change it for the number of gene tree that we define before
        ms_convert_res <- readLines("ms_convert_res.txt")
        pos_cut_ms_string <- regexpr(pattern="-T", text= ms_convert_res)[1]-3
        tail_ms_mod <- substring(ms_convert_res, pos_cut_ms_string + 2, nchar(ms_convert_res))
        head_ms_mod <- substring(ms_convert_res, 1, pos_cut_ms_string)
        
        # Generate gene tree with ms
        system(paste(head_ms_mod , k , tail_ms_mod, 
                     ">", filename,"_GT.txt", sep=""))
        
        net_counter <- net_counter + 1
        
        setwd("../")
        
        
        
      } #GAB hay un problema acá, no cierra el corchete
    }
  }
}


