### R script that simulates networks under the birth-death-hybridization
### model implemented in SiPhyNetworks using maximum granularity

### Original Parameters ----------------------------------------------------------

library(SiPhyNetwork) # library to simulate Networks, this use "ape" as dependence
# starting execution in scripts
source("00.generate_seeds.R")
source("functions.R") # load functions that operate on networks


# set the seed
set.seed(global_seed)

##make a seed file for the seeds used to generate the phylogeny
seed_file <- expand.grid(n_phy=1:n_phy,setting_no=setting_no, phy_seed= NA)

for(rw_no in setting_no){
  print(paste('simulating parameter setting',rw_no))
  par_setting <- pars[rw_no,]
  nu<- as.numeric(par_setting[2])
  ntips <- as.numeric(par_setting[3])
  level1 <- as.logical(par_setting[5])
  
  setting_folder <- paste(data_folder,"pars_",rw_no,'/',sep='')
  dir.create(setting_folder)
  
  n_success<- 1 #number of successful simulations for a given setting
  while(n_success <= n_phy){   #Keep going until we get enough valid networks
    seed<-sample(1e8,1)
    set.seed(seed)
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
    if (!is.phylo(network)) {
      next
    }
    # now check whether it is network or restart if not
    if (nrow(network$reticulation)==0) {
      next
    }
    
    ##TODO delete things above the lsa 
    
    
    ##check  2-cycles and 3-cycles
    good_net<-T
    for(hyb in network$reticulation[,2]){
      good_net <- check_2cycle(net,hyb) || check_3cycle(net,hyb)
      if(!good_net){ ##we found a  2 or 3 cycle
        break
      }
    }
    if(!good_net){
      next
    }
    
    
    # Find rows with both matching starting and end nodes
    is_match <- (from_match_indices == to_match_indices) & !is.na(from_match_indices)
    
    ##Make sure the network is level 1, if desired 
    if(!level1 && getNetworkLevel(network)==1){
      next
    }
    
    # if it is, write to file and generate/save seeds for downstream steps
    phy_folder <- paste(setting_folder,"net_",n_success,'/',sep='')
    dir.create(phy_folder)
    SiPhyNetwork::write.net(net = network, 
                            file = paste(phy_folder,"network.extnewick",sep=''))
    seed_file[((rw_no-1)*n_phy)+(n_success),3] <- seed
    
    ##create seeds for the downstream simulations/analyses
    sim_seeds <- expand.grid(rep = 1:n_reps,gt_seed=NA, snaq_seed = NA)
    sim_seeds$gt_seed   <- sample(1e8,size=nrow(sim_seeds))
    sim_seeds$snaq_seed <- sample(1e8,size=nrow(sim_seeds))
    write.csv(sim_seeds,file = paste(phy_folder,"seeds.csv",sep=''))
    
    n_success<- n_success+1
  }
  
}

##Generate CFs and starting tree
system(paste("julia ./02.sim_gts_calc_cf.jl",n_phy,n_reps))

##create folders for inference
job_no <- 0
dir.create('../jobs')
for(rw_no in setting_no){ #full compression
  setting_folder <- paste(data_folder,"pars_",rw_no,'/',sep='')
  for(phy_no in 1:n_phy){
    print(paste('rw',rw_no,'phy',phy_no))
    phy_folder <- paste(setting_folder,"net_",phy_no,'/',sep='')
    seed_file <- read.csv(paste(phy_folder,'seeds.csv',sep=''))
    for(rep_no in 1:n_reps){
      rep_folder <- paste(phy_folder,"rep_",rep_no,'/',sep='')
      
      ##generate starting tree -- NOW done in 02
      #command <- paste("tree-qmc --fast -i ",rep_folder,'gene_trees.newick -o ',
      #                 rep_folder, 'starting_tree.newick',sep='')
      #system(command,show.output.on.console = F)
      
      ##Generate data frame with the information for estimation
      est_pars<-data.frame(hmax = hmax,
                           snaq_seed = seed_file$snaq_seed[rep_no],
                           nthreads = nthreads,
                           nruns=nruns)
      write.csv(est_pars,paste(rep_folder,'est_pars.csv',sep=''))
      output_dir <- paste("out",sep='')
      dir.create(paste(rep_folder,output_dir,sep=''))
      
      ##Create folders with data to be sent to clusters.
      dir.create(paste('../jobs/job_',job_no,sep=''))
      job_files <-paste(paste(c('starting_tree.newick','CFs.csv','est_pars.csv'),sep=''),collapse = ' ')
      command <- paste("tar -czf ../../jobs/job_",job_no,'/files.tar.gz '," -C ",rep_folder," ",job_files," ",output_dir,sep='')
      system(command)
      unlink(paste(rep_folder,output_dir,sep='')) #don't keep this directory around
      job_no<- job_no+1
    }
  }
}


##link the job number to the parameter combo, phy_no, and rep_no
##Each row in job_map corresponds to the job numbers associated parameters
df_expanded <- do.call(rbind, lapply(1:nrow(pars), function(i) {
  cbind(pars[rep(i, n_phy), ], phy_no = 1:n_phy)
}))
df_expanded <- do.call(rbind, lapply(1:nrow(df_expanded), function(i) {
  cbind(df_expanded[rep(i, n_reps), ], rep_no = 1:n_reps)
}))
df_expanded$job_no <- 0:nrow(df_expanded)
write.csv(df_expanded,"../../data/job_map.csv",row.names = F)













