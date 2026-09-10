library(MSCquartets)
library(SiPhyNetwork)
library(ape)
library(TreeDist)
library(data.table)
library(foreach)
library(doParallel)
library(caret) 
library(rpart)
library(ipred)

# 1. Setup the parallel backend
# Use detectCores() - 1 to leave one core free for your operating system
num_cores <- detectCores() - 1
# outfile = "" ensures that your print() statements will still show up in the console
cl <- makeCluster(num_cores, outfile = "")
registerDoParallel(cl)

print(paste("Running on", num_cores, "cores..."))

# 2. Outer loop is back to a standard for-loop
for(par_no in 1:36) {
  
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
  
  ##create directories as needed 
  pars_dir <- paste("../summarized_results/pars_",par_no,"/",sep='')
  dir.create(pars_dir,recursive = T,showWarnings = F) # Changed to F to avoid clutter in parallel
  dir.create(paste(pars_dir,'tob',sep = ''),showWarnings = F,recursive=TRUE)
  
  tob_file <- paste(pars_dir,'tob.csv',sep='')
  missing_file <- paste(pars_dir,'missing.csv',sep='')
  tob_compat_file <- paste(pars_dir,'tob/compat.csv',sep='')
  
  #check if the summary files exist
  if(!file.exists(tob_file)){ #create empty csv if they do not
    write.csv(tob_dat,tob_file,row.names = F)
    write.csv(missing_dat,missing_file,row.names = F)
    write.csv(tob_compat_dat,tob_compat_file,row.names = F) 
  }
  
  ID_df <- read.csv(tob_file)
  done_keys <- paste(ID_df$phy, ID_df$rep, ID_df$hmax, sep = "_")
  
  # Parallelize the inner loop over phy_no instead
  phy_results <- foreach(phy_no = 1:150, .packages = c("MSCquartets", "SiPhyNetwork", "ape", "TreeDist", "data.table")) %dopar% {
    
    list_ind<-1
    tob_list <-list()
    missing_list<-list()
    compat_list <-list()
    
    print(paste("we are at par",par_no, "and phy",phy_no))
    #Compute things on the true network
    found_true_net <-FALSE
    for(try_net in 1:30){
      true_net_loc = paste("../output/pars/pars_",par_no,"/phy_",phy_no,"/rep_",try_net,"/network.extnewick",sep='')
      if(!file.exists(true_net_loc)){ #move on if the file does not exist
        next
      }else{
        found_true_net<-TRUE
        break
      }
    }
    if(!found_true_net){
      # Replaced 'next' with an empty return for the foreach loop
      return(list(tob = data.frame(), missing = data.frame(), compat = data.frame())) 
    }
    
    true_net <- read.net(true_net_loc)
    true_tob=treeOfBlobs(true_net,plot=F);
    d=degree(true_tob,details=T) 
    true_blob_nds  <- which(d>3)
    
    ##true Clades of blobs are:
    true_blobs<- list()
    for(blob_nd in true_blob_nds){
      true_blobs[[blob_nd]] <- nodeGroups(true_tob,blob_nd)
    }
    
    for(rep_no in 1:30){
      for(h_no in 1:5){
        
        current_key <- paste(phy_no, rep_no, h_no, sep = "_")
        
        if(current_key %in% done_keys) {
          next
        }
        
        #create empty dfs for each row
        rw_tob <- tob_dat
        rw_missing <- missing_dat
        rw_tob_compat <- tob_compat_dat
        
        est_net_loc = paste("../output/pars/pars_",par_no,"/phy_",phy_no,"/rep_",rep_no,"/h_",h_no,".out",sep='')
        
        if(!file.exists(est_net_loc)){ 
          rw_missing<-rbind(rw_missing,c(phy_no=phy_no,rep_no=rep_no,h_no=h_no))
          next
        }
        
        con <- file(est_net_loc,"r")
        first_line <- readLines(con,n=1)
        close(con)
        est_net <- read.net(text=first_line)
        if(!is.null(est_net$node.label)){
          empty_labels<-est_net$node.label==''
          est_net$node.label<-make.unique(est_net$node.label) 
          est_net$node.label[empty_labels]='' 
        }
        is_tree <- nrow(est_net$reticulation)==0
        
        if(!all(sort(est_net$tip.label) == sort(true_net$tip.label))){
          print(paste("there was a tip mismatch at",rep_no,sep=''))
          next
        }
        
        if(is_tree){
          est_net <- read.tree(text=first_line)
          est_tob <- est_net 
          
        }else{
          est_tob <-  treeOfBlobs(est_net,plot=F)
        }
        d=degree(est_tob,details=T) 
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
          for(true_blob_nd in true_blob_nds){
            true_blob_clades<-true_blobs[[true_blob_nd]]
            for(true_clade in true_blob_clades){
              found_clade <- F
              for(est_clade in est_blobs[[blob_nd]]){ 
                matching_tips <- true_tob$tip.label[true_clade] %in% est_tob$tip.label[est_clade]
                if(all(matching_tips)){ 
                  found_clade <- T
                  break  
                }
              }
              if(!found_clade){
                compat_blob <- F
                break
              }else{
                compat_blob <- T
              }
              if(!compat_blob){
                break 
              }
            }
            if(compat_blob){ 
              break 
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
        
        suppressMessages(rw_tob<-rbind(rw_tob,c(
          phy_no,rep_no,h_no,
          RF.dist(true_tob,est_tob),
          ClusteringInfoDist(true_tob,est_tob),
          MatchingSplitDistance(true_tob,est_tob),
          MatchingSplitInfoDistance(true_tob,est_tob)
        )))
        
        rw_tob_compat<-rbind(rw_tob_compat,compat_frame)
        
        colnames(rw_tob)<-c('phy','rep','hmax','rf_dist','clust_dist','split_dist','split_info_dist')
        rw_tob$phy <- as.integer(rw_tob$phy)
        rw_tob$rep <- as.integer(rw_tob$rep)
        rw_tob$hmax <- as.integer(rw_tob$hmax)
        colnames(rw_missing) <- c("phy",'rep','hmax')
        colnames(rw_tob_compat) <-c("phy",'rep','hmax','true_blob','est_blob')
        
        tob_list[[list_ind]]<-rw_tob
        missing_list[[list_ind]]<-rw_missing
        compat_list[[list_ind]]<-rw_tob_compat
        list_ind<-list_ind+1
        
      } #end hmax
    } #end rep_no
    
    phy_tob <- data.table::rbindlist(tob_list)
    phy_missing <-data.table::rbindlist(missing_list)
    phy_compat <- data.table::rbindlist(compat_list)
    
    # Return lists to the master process instead of writing inside parallel threads
    list(tob = phy_tob, missing = phy_missing, compat = phy_compat)
    
  } #end phy_no foreach loop
  
  # Safely combine results from all parallel workers for this par_no
  all_tob <- data.table::rbindlist(lapply(phy_results, function(x) x$tob))
  all_missing <- data.table::rbindlist(lapply(phy_results, function(x) x$missing))
  all_compat <- data.table::rbindlist(lapply(phy_results, function(x) x$compat))
  
  # Append to files safely, once per par_no
  if(nrow(all_tob) > 0) write.table(all_tob,tob_file,append = T, sep = ",",row.names = F,col.names = F)
  if(nrow(all_missing) > 0) write.table(all_missing,missing_file,append=T, sep = ",",row.names = F,col.names = F)
  if(nrow(all_compat) > 0) write.table(all_compat,tob_compat_file,append=T, sep = ",",row.names = F,col.names=F)
  
} #end par_no standard loop

# 3. Clean up the cluster
stopCluster(cl)
print("Finished parallel processing.")