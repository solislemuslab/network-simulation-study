###########################################
## Map job numbers to parameter settings ##
###########################################
if (!requireNamespace("foreach", quietly = TRUE)) install.packages("foreach")
if (!requireNamespace("doParallel", quietly = TRUE)) install.packages("doParallel")

library(foreach)
library(doParallel)


source("functions.R")



dir.create("../output/pars", recursive = TRUE, showWarnings = FALSE)
dir.create("../output/all", recursive = TRUE, showWarnings = FALSE)


if(file.exists("../output.tar.gz")){
  print("Extracting main output archive...")
  system2("tar", args = c("-xzf", "../output.tar.gz", "-C", "../output/all"))
  dir.create("../processed_outputs", recursive = TRUE, showWarnings = FALSE)
  
  
  timestamp <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
  new_filename <- paste0("../processed_outputs/output_", timestamp, ".tar.gz")
  file.rename("../output.tar.gz", new_filename)
  print(paste("Archive moved and renamed to:", new_filename))
}


num_cores <- max(1, parallel::detectCores() - 1) 
cl <- makeCluster(num_cores,outfile='')
registerDoParallel(cl)

total_jobs <- 30 * 150 * 36

results <- foreach(i = 0:total_jobs, .combine = rbind, .packages = c()) %dopar% {
  
  source("functions.R")
  pars <-ID2pars(i)
  setting_no <-pars[1]
  phy_no <- pars[2]
  rep_no <- pars[3]

  # Print a progress heartbeat every 1000 jobs
  if (i %% 1000 == 0) {
    cat(sprintf("Worker processing job %d...\n", i))
  }
  
  job_dir <- file.path("../output/all/jobs", paste0("job_", i))
  emp_dir <- file.path(job_dir, paste0("emp_est_out", i, ".tar.gz"))
  files_file <- file.path(job_dir, "files.tar.gz")
  
  pars_dir <- file.path("../output/pars", paste0("pars_", setting_no), paste0("phy_", phy_no))
  rep_dir <- file.path(pars_dir, paste0("rep_", rep_no))
  
  
  
  # Skip logic
  if (dir.exists(rep_dir)) {
    return(data.frame(id = i, status = "happy"))
  }
  
  if (!file.exists(emp_dir)) {
    return(data.frame(id = i, status = "skipped"))
  }
  
  # Create directory
  dir.create(rep_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Use list.files ONCE
  search_pattern <- paste0(".*", "\\.tar\\.gz$")
  all_tarfiles <- list.files(path = job_dir, pattern = search_pattern, full.names = TRUE)
  
  # Find and remove CF_out files
  cf_files <- grep("/CF_out", all_tarfiles, value = TRUE)
  if (length(cf_files) > 0) {
    file.remove(cf_files)
    # Filter the list in RAM rather than calling list.files() again
    all_tarfiles <- setdiff(all_tarfiles, cf_files) 
  }
  
  
  # Decompress using system2 (bypasses cmd.exe overhead) and NO verbose '-v' flag
  for (tarfile in all_tarfiles) {
    system2("tar", args = c("-xzf", tarfile, "-C", rep_dir), stdout = FALSE, stderr = FALSE)
  }
  
  #if (file.exists(files_file)) {
  #  system2("tar", args = c("-xzf", files_file, "-C", rep_dir), stdout = FALSE, stderr = FALSE)
  #}
  
  return(data.frame(id = i, status = "processed"))
}
# Stop the parallel cluster
stopCluster(cl)

