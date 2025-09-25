### R script that simulates networks under the birth-death-hybridization
### model implemented in SiPhyNetworks using maximum granularity

##using the networks generated from data par1

data_folder <- "../data_toy/"
##Generate CFs and starting tree
system(paste("julia ./02.sim_gts_toy.jl"))

##Generate starting trees and create folders for inference
job_no <- 0
dir.create('../jobs_toy')

phy_nos <- c(1,5,13,15,16,18,25,34)
n_reps <- 30
set.seed(48284)
for(phy_no in phy_nos){
  phy_folder <- paste(data_folder,"net_",phy_no,'/',sep='')
  for(rep_no in 1:n_reps){
    rep_folder <- paste(phy_folder,"rep_",rep_no,'/',sep='')
    

    
    ##Generate data frame with the information for estimation
    est_pars<-data.frame(hmax = 5,
                         snaq_seed = sample(1e8,1),
                         nthreads = 10,
                         nruns=5)
    write.csv(est_pars,paste(rep_folder,'est_pars.csv',sep=''))
    output_dir <- paste("out",sep='')
    dir.create(paste(rep_folder,output_dir,sep=''))
    
    ##Create folders with data to be sent to clusters.
    dir.create(paste('../jobs_toy/job_',job_no,sep=''))
    job_files <-paste(paste(c('starting_tree.newick','CFs.csv','est_pars.csv'),sep=''),collapse = ' ')
    command <- paste("tar -czf ../jobs_toy/job_",job_no,'/files.tar.gz '," -C ",rep_folder," ",job_files," ",output_dir,sep='')
    system(command)
    unlink(paste(rep_folder,output_dir,sep='')) #don't keep this directory around
    job_no<- job_no+1
  }
}

