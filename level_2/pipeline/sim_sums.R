library(ape)
library(SiPhyNetwork)
library(phangorn)
library(MSCquartets)
source("functions.R")
source("00.generate_seeds.R") #get pars,n_phy

data_dir <- "../data/"

sum_dat <- matrix(NA,nrow = nrow(pars)*n_phy,ncol = 12 )

sum_dat <- data.frame(sum_dat)
colnames(sum_dat) <- c("par","phy","nrets","level","blobs","nt_blobs","TC","TB","FU","cycle2","cycle3","LSA")

total_row<-1
for(i in 1:nrow(pars)){
  rw<-pars[i,]
  print(paste("setting number",rw$setting_no))
  rw<-pars[i,]
  par_dir <- paste(data_dir,
                   'pars_',rw$setting_no,'/',
                   sep='')

  
  for(phy_no in 1:n_phy){
    net <- read.net(paste(par_dir,"net_",phy_no,"/network.extnewick",sep=''))
    
    ##Compute the biconnected components for the level and number of blobs
    rt <- as.integer(length(net$tip.labels) + 1)
    edges <- rbind(net$edge, net$reticulation)
    hyb_nds <- net$reticulation[, 2]
    mode(edges) <- "integer"
    nNode <- length(net$tip.label) + net$Nnode
    blobs <- biconnectedComponents(edges, rt, nNode)
    net_level <- 1
    for (blob in blobs) {
      blob <- blob + 1
      blob_nds <- unique(as.vector(edges[blob, ]))
      blob_nds <- blob_nds[blob_nds %in% hyb_nds]
      net_level <- max(c(net_level, length(blob_nds)))
    }
    
    true_tob=treeOfBlobs(net,plot=F);
    d=degree(true_tob,details=T) # This function produces a table telling how many nodes of degree n  we want  the number of nodes that have degree greater than 3
    blob_nds  <- which(d>3)
    blobs_nontrivial <-length(blob_nds)
    
    ret_nodes <- net$reticulation[,2]
    
  
    # ##get the branch lengths
    # net$node.label<-NULL
    # 
    # internal_e <- net$edge.length[!(net$edge[,2] %in% 1:length(net$tip.label))]
    # bt<- branching.times(net)
    # hyb_e = bt[as.character(net$reticulation[,1])]-bt[as.character(net$reticulation[,2])]
    # e <- c(internal_e,hyb_e)
    # e <- e[1e-9<e] ## get rid of branches with effectively zero length. rounding error
    
    ##Get the 2cycle,3cycle, and non-LSA networks
    ret_nds<- unique(net$reticulation[,2])
    edges <- rbind(net$edge,net$reticulation)
    
    #look for 2 and 3 cycles
    two_cycles<-0
    three_cycles<-0
    for(i in 1:nrow(net$reticulation)){
      nd <- net$reticulation[i,2]
      par_nds <- edges[which(edges[,2]==nd),1]
      if(par_nds[1]==par_nds[2]){ #both edges have the same parent
        #we have a 2cycle
        two_cycles <- two_cycles+1
        next
      }
      
      ##look for 3cycle
      par1_par <- edges[which(edges[,2]==par_nds[1]),1] #parent of the left hybrid parent
      par2_par <- edges[which(edges[,2]==par_nds[2]),1] #parent of the right hybrid parent
      
      #three cycle if
      #the parent of the left hybrid parent is the same as the right hybrid parent
      #OR
      #the parent of the right hybrid parent is the same as the left hybrid parent
      if(all(par1_par==par_nds[2]) || all(par2_par==par_nds[1])){
        three_cycles<- three_cycles+1
      }
      
    }
    
    
    #look for LSA after root
    #ltt <-ltt.network(net,bt)
    #max_rw <-1:(nrow(ltt)) < which(ltt$n_lineages==length(net$tip.label))[1] #only works for SSA(?)
    #ltt <- ltt[max_rw,] # only keep proper rows at the end
    
    #LSA<-any(ltt$n_lineages==1)
    
    ##save the 
    rw_dat<-  c(rw$setting_no,
                phy_no,
                nrow(net$reticulation),
                net_level,length(blobs),blobs_nontrivial,
                isTreeChild(net),
                isTreeBased(net),
                isFUstable(net),
                two_cycles,
                three_cycles,
                0)
    sum_dat[total_row,]<-rw_dat
    total_row<-total_row+1
  }
}
write.csv(sum_dat,'../data/sim_net_props.csv',row.names = F)


