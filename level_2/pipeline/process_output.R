###########################################
## Map job numbers to parameter settings ##
###########################################
source("functions.R")
skipped<-c() ## 
happy <- c()
dir.create("../output")
dir.create("../output/pars")
dir.create("../output/all")
shell("tar -xvzf ../output.tar.gz -C ../output/all")

for(i in 0:(30*150*36)){
  pars <-ID2pars(i)
  setting_no <-pars[1]
  phy_no <- pars[2]
  rep_no <- pars[3]
  if((i %% 100)==0) {
    print(i)
  }
  job_dirs <- paste("../output/all/jobs/job_",i,"/*_out",i,".tar.gz",sep='')
  emp_dir <-paste("../output/all/jobs/job_",i,"/emp_est_out",i,".tar.gz",sep='')
  files_file <- paste("../output/all/jobs/job_",i,"/files",".tar.gz",sep='')
  
  pars_dir <- paste("../output/pars/pars_",setting_no,
                    "/phy_",phy_no,
                    sep='')
  dir.create(pars_dir,recursive = T,showWarnings = F)
  
  #don't do things if we already pasted or
  #if the output doesn't exist
  rep_dir <- paste(pars_dir,
                   "/rep_",rep_no,
                   sep='')
  if(dir.exists(rep_dir)){#we already moved things
    happy <- c(happy,i)
    #print(paste('job',i,'already was transferred'))
    next
  }
  if(!file.exists(emp_dir)){ #output doesn't exit
    #print(paste('job',i,'was skipped'))
    skipped<-c(skipped,i) #record what we skipped
    next
  } 
  

  
  # Extract the directory (dirname) and the pattern (basename)
  target_dir <- dirname(job_dirs)
  search_pattern <- basename(job_dirs)
  
  # Use list.files to find all matching files
  all_tarfiles <- list.files(
    path = target_dir,
    pattern = glob2rx(search_pattern), # Use glob2rx to convert shell glob to regex
    full.names = TRUE # Crucial: gives the full path to each file
  )
  
  target_file <- grep("/CF_out", all_tarfiles, value = TRUE)
  if(length(target_file)==0){
    #next #move on, nothing to see here
  }else {
    print(paste('we destroy CF in ',i,sep=''))
    #there are CF tar files. destroy them
    for(target in target_file){
      #destroy
      file.remove(target)
    }
  }
  # Use list.files to find all matching files
  all_tarfiles <- list.files(
    path = target_dir,
    pattern = glob2rx(search_pattern), # Use glob2rx to convert shell glob to regex
    full.names = TRUE # Crucial: gives the full path to each file
  )
  
  #re-extract everything
  dir.create(rep_dir)
  for(tarfile in all_tarfiles){
    shell(paste(
      "tar -xvzf",tarfile,"-C",
      rep_dir,sep= ' '),ignore.stdout=T)
  }
  shell(paste(
    "tar -xvzf",files_file,"-C",
    rep_dir,sep= ' '),ignore.stdout=T)

}
write.csv(data.frame(skipped),"missing_jobs.csv",row.names = F)
