### R script that simulates networks under the birth-death-hybridization
### model implemented in SiPhyNetworks using maximum granularity

### Original Parameters ----------------------------------------------------------

library(SiPhyNetwork) # library to simulate Networks, this use "ape" as dependence
library(ape)
# starting execution in scripts
source("00.generate_seeds.R")
source("functions.R") # load functions that operate on networks


# set the seed


set.seed(global_seed)
setting_seeds <-sample(1e8,length(setting_no))
##make a seed file for the seeds used to generate the phylogeny
seed_file <- expand.grid(n_phy=1:n_phy,setting_no=setting_no, phy_seed= NA)

dir.create('../jobs',showWarnings = F)

for(rw_no in setting_no[1:36]){
  set.seed(setting_seeds[rw_no])
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
                                hyb.inher.fxn = make.beta.draw(10, 10),
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
    n_rets<- nrow(network$reticulation)
    # now check whether it is network or restart if not
    if (n_rets==0) {
      next
    }
    
    ##check  2-cycles and 3-cycles
    has_bad_cycle <- check_2cycle(network) || check_3cycle(network)
    if(has_bad_cycle){
      next
    }
    
    #make sure the root is the LSA 
    ltt<-ltt.network(network)
    if(any(ltt$n_lineages ==1)){ #we only get back down to one lineage after a hybridization. Root not LSA
      next
    }

    ##Make sure the network is level 1, if desired 
    if(level1){
      if(getNetworkLevel(network)!=1){
        next
      }
    }else{
      if(getNetworkLevel(network)==1){
        next
      }
    }

    
    # if it is, write to file and generate/save seeds for downstream steps
    phy_folder <- paste(setting_folder,"net_",n_success,'/',sep='')
    dir.create(phy_folder)
    SiPhyNetwork::write.net(net = network, 
                            file = paste(phy_folder,"network.extnewick",sep=''))
    
    
    
    seed_file[((rw_no-1)*n_phy)+(n_success),3] <- seed
    print( paste("parameter", rw_no,"and phy",n_success))
    for(rep_no in 1:n_reps){#generate a jobs folder for each rep
      rep_folder <- paste(phy_folder,"rep_",rep_no,'/',sep='')
      dir.create(rep_folder,showWarnings = F)
      
      SiPhyNetwork::write.net(net = network, 
                              file = paste(rep_folder,"network.extnewick",sep=''))
      ##create seeds for the downstream simulations/analyses
      est_pars<-data.frame(hmax = min(c(n_rets,5)),
                       snaq_seed = sample(1e8,1),
                       gt_seed = sample(1e8,1),
                       gof_seed = sample(1e8,1),
                       nruns=nruns,
                       ngt= par_setting$ngt,
                       par_no = par_setting,
                       phy = n_success,
                       rep = rep_no
                       )
      
      write.csv(est_pars,file = paste(rep_folder,"seed_pars.csv",sep=''))
      
      
      ##create the compressed job folder
      job_no <- format(get_jobID(rw_no,n_success,rep_no),scientific=F)
      output_dir <- paste("out",sep='')
      dir.create(paste(rep_folder,output_dir,sep=''))
      
      ##Create folders with data to be sent to clusters.
      dir.create(paste('../jobs/job_',job_no,sep=''))

      job_files <-paste(paste(c('network.extnewick','seed_pars.csv'),sep=''),collapse = ' ')
      command <- paste("tar -czf ../jobs/job_",format(job_no,scientific=F),'/files.tar.gz '," -C ",rep_folder," ",job_files," ",output_dir,sep='')
      system(command)
      unlink(paste(rep_folder,output_dir,sep='')) #don't keep this directory around
      

      
    }

    n_success<- n_success+1
  }
}






