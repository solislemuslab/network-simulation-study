###########################################
## Map job numbers to parameter settings ##
###########################################
job_map <-read.csv("../data/job_map.csv")
skipped<-c() ## 
happy <- c()
for(i in 1:nrow(job_map)){
  rw <-job_map[i,]
  if((i %% 500)==0) {
    print(i)
  }
  job_dir <- paste("../output/all_outputs/output_",rw$job_no,"/out",sep='')

  
  pars_dir <- paste("../output/pars/pars_",rw$setting_no,
                    "/phy_",rw$phy_no,
                    sep='')
  dir.create(pars_dir,recursive = T,showWarnings = F)
  
  #don't do things if we already pasted or
  #if the output doesn't exist
  rep_dir <- paste(pars_dir,
                   "/rep_",rw$rep_no,
                   sep='')
  if(dir.exists(rep_dir)){#we already moved things
    happy <- c(happy,rw$job_no)
    print(paste('job',rw$job_no,'already was transferred'))
    next
  }
  if(!dir.exists(job_dir)){ #output doesn't exit
    print(paste('job',rw$job_no,'was skipped'))
    skipped<-c(skipped,rw$job_no) #record what we skipped
    next
  }
  file.rename(from=job_dir,to=rep_dir)
}
write.csv(data.frame(skipped),"missing_jobs.csv",row.names = F)
