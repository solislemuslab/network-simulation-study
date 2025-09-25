library(MSCquartets)
library(SiPhyNetwork)
library(ape)
library(TreeDist)

for(par_no in 1:36){
  tob_dat <- data.frame(phy = integer(),
                        rep = integer(),
                        hmax = integer(),
                        rf_dist = numeric(),
                        clust_dist = numeric(),
                        split_dist = numeric(),
                        split_info_dist =numeric()
  )
  missing_dat <- data.frame(phy_no = integer(),
                            rep_no = integer(),
                            hmax = integer()
    
  )
  tob_compat_dat <-data.frame(phy=integer(),
                              rep=integer(),
                              hmax=integer(),
                              true_blob=integer(),
                              est_blob=integer())
  
  dir.create(paste("../summarized_results/pars_",par_no,'/tob',sep = ''),showWarnings = F)
  
for(phy_no in 1:150){
  #Compute things on the true network
  true_net_loc = paste("../data/pars_",par_no,"/net_",phy_no,"/network.extnewick",sep='')
  true_net <- read.net(true_net_loc)
  true_tob=treeOfBlobs(true_net,plot=F);
  d=degree(true_tob,details=T) # This function produces a table telling how many nodes of degree n  we want  the number of nodes that have degree greater than 3
  true_blob_nds  <- which(d>3)
  ##true Clades of blobs are:
  true_blobs<- list()
  for(blob_nd in true_blob_nds){
    true_blobs[[blob_nd]] <- nodeGroups(true_tob,blob_nd)
  }
  
  print(paste("we are at par",par_no, "and phy",phy_no))
  for(rep_no in 1:30){
  for(h_no in 1:5){
    
    est_net_loc = paste("../output/pars/pars_",par_no,"/phy_",phy_no,"/rep_",rep_no,"/h_",h_no,".out",sep='')
    
    if(!file.exists(est_net_loc)){ #move on if the file does not exist
      missing_dat<-rbind(missing_dat,c(phy_no=phy_no,rep_no=rep_no,h_no=h_no))
      next
    }
    
    con <- file(est_net_loc,"r")
    first_line <- readLines(con,n=1)
    close(con)
    est_net <- read.net(text=first_line)
  
    est_tob <- treeOfBlobs(est_net,plot=F)
    d=degree(est_tob,details=T) # This function produces a table telling how many nodes of degree n  we want  the number of nodes that have degree greater than 3
    est_blob_nds  <- which(d>3)
    
    
    ##est clades of blobs are:
    est_blobs <- list()
    compatible_blobs <- vector("list", max(c(est_blob_nds,0)))
    
    compat_frame<-data.frame(phy=integer(),rep=integer(),hmax=integer(),
      true_blob=integer(),est_blob=integer())
    for(blob_nd in est_blob_nds){
      est_blobs[[blob_nd]] <- nodeGroups(est_tob,blob_nd)
      
      #check compatibility of est blob clades with the true one
      compat_blob<-F
      for(true_blob_nd in true_blob_nds){##Look over all true blob clades
        ##Check compatibility between true_blob_clade and est_blob_clades[[blob_nd]]
        ##Compatibility means that any clades defined in the true are:
        #at least found together, but can be lumped with other clades
        true_blob_clades<-true_blobs[[true_blob_nd]]
        #print(paste('we are looking for node:',true_blob_nd))
        for(true_clade in true_blob_clades){
          ##Check if this clade is found wholly somewhere among the est_clades
          found_clade <- F
          #print(paste('we are looking for clade',paste(true_tob$tip.label[true_clade],collapse='_')))
          for(est_clade in est_blobs[[blob_nd]]){ #look at each est clade
            matching_tips <- true_tob$tip.label[true_clade] %in% est_tob$tip.label[est_clade]
            #print(true_tob$tip.label[true_clade])
            #print(est_tob$tip.label[est_clade])
            if(all(matching_tips)){ 
              ## This clade was found!
              found_clade <- T
              break ## stop looking for this clade 
            }
          }
          if(!found_clade){## we couldn't find the clade
            #this true blob is not compatible 
            compat_blob <- F
            
            #print(paste('we couldnt find clade',paste(est_tob$tip.label[est_clade],collapse=' '),' and are in blob node:',blob_nd))
            break
          }else{
            compat_blob <- T
          }
          if(!compat_blob){
            break #the blob is not compatible, no need to keep looking at true clades
          }
        }
        if(compat_blob){ ##if we made it this far and it is compatible
         break ## we can stop looking for compatibility and move on
        }
      }
      if(compat_blob){
        compatible_blobs[[blob_nd]]<- c(compatible_blobs[[blob_nd]],true_blob_nd)
        new_rw <-data.frame(phy=phy_no,rep=rep_no,hmax=h_no,
                            true_blob=true_blob_nd,est_blob=blob_nd)
        compat_frame<-rbind(compat_frame,new_rw)
      }else{
        new_rw <-data.frame(phy=phy_no,rep=rep_no,hmax=h_no,
                            true_blob=-1,est_blob=blob_nd)
        compat_frame<-rbind(compat_frame,new_rw)
      }
    }
    
    
    
    
    suppressMessages(tob_dat<-rbind(tob_dat,c(
                   phy_no,rep_no,h_no,
                   RF.dist(true_tob,est_tob),
                   ClusteringInfoDist(true_tob,est_tob),
                   MatchingSplitDistance(true_tob,est_tob),
                   MatchingSplitInfoDistance(true_tob,est_tob)
                   )))
    # write.csv(compat_frame,paste("../summarized_results/pars_",par_no,'/compat_blobs/phy_',phy_no,'_rep_',rep_no,'_h_',h_no,'.csv',sep=''),row.names = F)
    tob_compat_dat <- rbind(tob_compat_dat,compat_frame)
  } #end hmax
} #end rep_no
} #end phy_no
  colnames(tob_dat)<-c('phy','rep','hmax','rf_dist','clust_dist','split_dist','split_info_dist')
  tob_dat$phy <- as.integer(tob_dat$phy)
  tob_dat$rep <- as.integer(tob_dat$rep)
  tob_dat$hmax <- as.integer(tob_dat$hmax)
  
  colnames(missing_dat) <- c("phy",'rep','hmax')
  colnames(tob_compat_dat) <-c("phy",'rep','hmax','true_blob','est_blob')
  
  dir.create(paste("../summarized_results/pars_",par_no,sep = ''),showWarnings = F)
  write.csv(tob_dat,paste("../summarized_results/pars_",par_no,'/tob.csv',sep=''),row.names = F)
  write.csv(missing_dat,paste("../summarized_results/pars_",par_no,'/missing.csv',sep=''),row.names = F)
  write.csv(tob_compat_dat,paste("../summarized_results/pars_",par_no,'/tob/compat.csv',sep=''),row.names = F) 
} #end par_no




