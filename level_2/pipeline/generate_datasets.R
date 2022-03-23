


##################################################################GAB Parameters were specified on notion in https://www.notion.so/Project-on-SNaQ-limitations-80e419d6d58a42308e4bf0f5cceb9752#e5df0d1a4c97433ca14f2e6388be966a
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



script_path <-"/home/carlos/Desktop/Borrar/Experiment1"

setwd(script_path)


library(SiPhyNetwork)#library to simulate Networks, this use "ape" as dependence
#GAB There is something really odd with geiger::is.phylo: It is not intended for external user or indocumented (see ?is.phylo)

#GAB if we do not have alternative, we should 1) use class(tree) which should return "phylo" for such case, or try to copy original geiger's code. There's no point in adding a dependency only because of a function which is not documented or intended for users

#GAB See geiger's original code (https://github.com/mwpennell/geiger-v2/blob/6dffb87328f5806e02e7ce06d2b2d409b6182db2/R/utilities-phylo.R#L460). All what is.phylo does is "phylo" %in% class(x) where x is a tree.
# library(geiger)#is.phylo()##Information about tree goes extinct=0 and no extinct tips are sampled=1

#GAB create a directory for storing the networks and their subdirs. the .. is necessary as it assumes that the script generate_datasets.R us run from the directory pipeline, which should only contain code and not data
dir.create(path = "../data")


#3.1 Simulate networks
#GAB set each seed as below, They are all pseudo-random, but dependent on the first set.seed(2022) above

##General seed
set.seed(500)

#GAB nested for for generating the params one at a time for those which vary, i.e., ntips, nu, ngt

for (i in ntips) {
    for (j in nu) {
        for (k in ngt) {
			
			setwd("../data")
            # networks_raw will store the 150 networks for each combination of params
			r_seed <- sample.int(n = 1e6, size = 1)
			set.seed(r_seed)
	    #i=30		
	    #j=nu[1]
	    #k=ngt[1]		
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
            #GAB code for picking only the networks of interest, Networks-only, Ultrametric-only, We are now interested in taking only the good nets onwards
            ######## CODE HERE #######
			
			file_networks <- c()
			for(x in 1:length(networks)){
			  neet_i <- networks[[x]]
			  ret_numb <- nrow(neet_i$reticulation)
			  if (ret_numb>=1){
			  file_networks <- c(file_networks,x)
			  }
			}
			networks1<-networks[file_networks]
			
			
			################################
			#need to improve this part because when the string is so long it can not pass to julia
			#one alternative is split in suitable parts and pas to julia
			###############################			
			rnet_for_julia0<-write.net(networks1)
			
			rnet_for_julia<-paste0(rnet_for_julia0, collapse = "_")
			setwd(script_path)
			
			net_hyla_format0<-system(paste("extnewick2hybridlambda.jl ", "'", rnet_for_julia, "'", sep = ""),intern = T)

			net_hyla_format <- strsplit(net_hyla_format0, split = "_")[[1]]
			###################################################################
			
			
			
			
			setwd("../data")
			
            #GAB after selecting only the networks that we need, run the julia script extnewick2hybridlambda.jl for format conversion over each element network in networks
            net_counter <- 1
            for (l in 1:length(networks1)) {
				
				network<-networks1[l]
  
				#kl=net_hyla_format[l]
				hyla_seed = sample.int(n = 1e6, size = 1)
				
				gt_filename <- paste("genetree", net_counter, "_ntips", i, "_nu", j, "_ngt", k, "_rseed", r_seed, "_gtseed",hyla_seed,sep="")
								
				w<-system(paste("hybrid-Lambda -spcu ", "'",net_hyla_format[l],"'", " -num ", k," -seed ", hyla_seed, " -o ",gt_filename," > ult_test.txt 2>&1",sep=""))

				ultrametric_logic <- grep(x=readLines("ult_test.txt")[3],pattern="Dot figure generated in file")
				net_counter <- net_counter + 1

				if(length(ultrametric_logic)==0){
				  #GAB the network file in extnewick will be called appending the parameter values as well as a counter for numbering each network from 1 no length(networks)
				  network_filename <- paste("net", net_counter, "_ntips", i, "_nu", j, "_ngt", k, "_rseed-", r_seed)
				  # create the subdirectory for a single network using network_filename
				  dir.create(network_filename)
				  
				  # write the network to a file in extnewick format, inside the directory network_filename
				  rnet_name <- paste(network_filename, "/", network_filename, ".extnewick", sep ="")
				  SiPhyNetwork::write.net(net = network, file = rnet_name)

				  file.remove("ult_test.txt")
				  file.copy(paste(gt_filename,"_coal_unit",sep=""),network_filename)
				  file.remove(paste(gt_filename,"_coal_unit",sep=""))
				}else{
					  file.remove("ult_test.txt")
					  file.remove(paste(gt_filename,"_coal_unit",sep=""))
					}
				
				}
			
}}}
