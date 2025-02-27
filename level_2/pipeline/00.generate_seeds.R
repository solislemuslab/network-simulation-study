data_folder<- "../data/"
dir.create(data_folder)

##Constant parameters
lambda <- 0.9              ## speciation rate
mu <- 0                    ## extinction rate
hybprops <- c(1, 1, 1)     ## probabilities for each type of hybridization THIS NEEDS TO BE PROVIDED AS A STRING SEPARATED BY COMMAS
global_seed <- 2022         ## global seed


##simulation replication parameters
n_phy  <- 150 # number of phylogenies that we simulate for each setting
n_reps <- 30 #number of times we simulate gene trees for a given phylogeny and simulation setting

##Estimation parameters
hmax <- 5
nthreads <- 10
nruns <- 10


##Experiment Parameters
nu <- c(0.02, 0.04)        ## hybridization rate
ntips <- c(15,20, 30)     ## number of leaves in the network
ngt <- c(100, 1000, 10000) ## number of gene trees to simulate per network
level1 <-c(0,1)

pars<-expand.grid(nu=nu,ntips=ntips,ngt=ngt,level1=level1)
setting_no <- 1:nrow(pars)
pars<-cbind(setting_no,pars)

write.csv(pars,file = paste(data_folder,"pars.csv",sep=''))

