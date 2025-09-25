library(tidyverse)
library(ggplot2)
library(survival)
source("./functions.R")

setwd("../pipeline")

net_dat <- read.csv('../data/sim_net_props.csv')

par_no<-1

sim_folders<-c('branch_lengths','rets','level','blobs')
simfig_dir <- paste('../summarized_results/pars_',par_no,'/sims/figs/',sep = '')
dir.create(simfig_dir,recursive = T)


#######################################
#### Generate simulation summaries ####
#######################################

bl <- read.csv(paste('../data/pars_',par_no,'/branch_lengths.csv',sep=''))
prob_not_coal <- round(mean(exp(-bl$bl)),digits = 3)
hist_plot <- ggplot(bl, mapping = aes(x=bl))+geom_histogram()+labs(x='branch length',title = paste('Average probability of not coalescing: ',prob_not_coal,sep=''))
ggsave(filename=paste(simfig_dir,"hist.png",sep=''),plot=hist_plot,device = 'png')

par_dat<- net_dat %>% filter(par==par_no)

level_plot<-ggplot(par_dat,mapping = aes(x=level))+geom_histogram()+labs("Network Level")
ggsave(filename=paste(simfig_dir,"level.png",sep=''),plot=level_plot,device = 'png')

blob_plot <- ggplot(par_dat,mapping = aes(x=blobs))+geom_histogram()+labs("Number of blobs")
ggsave(filename=paste(simfig_dir,"blobs.png",sep=''),plot=blob_plot,device = 'png')

ret_plot <- ggplot(par_dat,mapping = aes(x=nrets))+geom_histogram()+labs("Number of reticulations")
ggsave(filename=paste(simfig_dir,"rets.png",sep=''),plot=ret_plot,device = 'png')


sum(par_dat$nt_blobs!=par_dat$blobs)


#########################################
#########################################
#####                              ######
#####       Summarize results      ######
#####                              ######
#########################################
#########################################


res_file_loc <- paste("../summarized_results/pars_",par_no,"/",sep='')
cluster_res <- read.csv(paste(res_file_loc,'clusters.csv',sep=''))
tob_res     <-read.csv(paste(res_file_loc,'tob.csv',sep=''))
res<- merge(cluster_res,tob_res,by=c('phy','rep','hmax'))

par_dat <- net_dat %>% ##data from the true network
  filter(par==par_no) %>%
  dplyr::select(phy,nrets,level,blobs,nt_blobs) 

res <- merge(res,par_dat,by='phy',all.x=TRUE) #add true phy info

hyb_res <- read.csv(paste("../summarized_results/pars_",par_no,"/found_clades.csv",sep=''))
hyb_res <- merge(hyb_res,par_dat,by='phy',all.x=TRUE) #add true phy info

tob_res <- read.csv(paste("../summarized_results/pars_",par_no,"/tob/compat.csv",sep=''))
tob_res <- merge(tob_res,par_dat,by='phy',all.x=TRUE) #add true phy info

dat_filters = list(
  list(name='hmax1',filter = expr(filter(hmax==1))),
  list(name='hmax2',filter = expr(filter(hmax==2))),
  list(name='hmax3',filter = expr(filter(hmax==3))),
  list(name='hmax4',filter = expr(filter(hmax==4))),
  list(name='hmax5',filter = expr(filter(hmax==5))),
  list(name='ntblob',filter = expr(filter(min(nt_blobs,5) == hmax))),
  list(name='blob',filter = expr(filter(min(blobs,5) == hmax))),
  list(name='level',filter = expr(filter(min(level,5) == hmax))),
  list(name='nrets',filter = expr(filter(min(nrets,5) == hmax)))
)
filter_names <- sapply(dat_filters, function(x) x$name)

hyb_dect_names<-c("dist50","externaltrue","blob_level","cycle_size")
hyb_pval_names<-paste('hyb_detection/',hyb_dect_names,'/pval',sep='')
hyb_coef_names<-paste('hyb_detection/',hyb_dect_names,'/coef',sep='')

plot_graphs<-c('PPV','TPR','f1',
               'true_clus_broad','true_clus_exact',
               'est_clus_broad','est_clus_exact',
               'rf',
               hyb_pval_names,
               hyb_coef_names,
               'tob_compat'
               )
plot_dirs<- paste('../summarized_results/pars_',par_no,'/figs/',plot_graphs,'/',sep='')
sapply(plot_dirs, dir.create, recursive = TRUE, showWarnings = FALSE)

# Create an empty data frame with those row names
col_names <-c(
  "hwcd",
  "broad_compat","exact_compat",
  "rf_dist","clust_dist","split_dist","split_info_dist",
  "ppv","tpr","f1",
  "broad_coverage",  "narrow_coverage","exact_coverage",
  "broad_mapping","exact_mapping",  
  "perc_compat"    
)
across_filter <- data.frame(matrix(ncol = length(col_names), nrow = 0))
for( i in seq_along(dat_filters)){
  print(i)
  dat_filter <- dat_filters[[i]]
  file_name  <- dat_filter$name
  filt       <- dat_filter$filter
  
  print(file_name)
  #filter the data based on 'filt'
  dat <-     eval(expr(res %>% group_by(phy) %>%
                         !!filt %>% ungroup()))
  hyb_dat <- eval(expr(hyb_res %>% group_by(phy) %>%
                         !!filt %>% ungroup()))
  tob_dat <- eval(expr(tob_res %>% group_by(phy) %>%
                         !!filt %>% ungroup()))
  
  
 
  suppressWarnings(suppressMessages(source('filtered_sum.R')))
  across_filter<-rbind(across_filter,cbind(across_phy,across_tob))
}


across_filter$name <- filter_names

hmax_dat<-across_filter[1:5,]



#############################
### Tree of Blobs Summary ###
#############################

par_no<-1

res_file_loc<- paste("../summarized_results/pars_",par_no,"/tob.csv",sep='')
res <-read.csv(res_file_loc)

summed_res <- res %>% group_by(hmax) %>%
  summarise(
            prop_exact_tob = sum(rf_dist == 0)/n(),
            rf_dist = mean(rf_dist),
            clust_dist = mean(clust_dist),
            split_dist = mean(split_dist),
            split_info_dist = mean(split_info_dist)
            )

blob_res <- res %>% group_by(phy) %>%
  filter(min(n_blobs,5) == hmax) %>% ungroup()

blob_res <- blob_res %>%
  summarise(
    prop_exact_tob = sum(rf_dist == 0)/n(),
    rf_dist = mean(rf_dist),
    clust_dist = mean(clust_dist),
    split_dist = mean(split_dist),
    split_info_dist = mean(split_info_dist)
  )

###################################
##### Goodness of Fit Summary #####
###################################


par_no<-13

res_file_loc<- paste("../summarized_results/pars_",par_no,"/GoF.csv",sep='')
res <-read.csv(res_file_loc)

summed_res <- res %>% group_by(hmax) %>%
  summarise(
    ave_pval = mean(pval),
    var_pval = var(pval),
    pval_2_5 = quantile(pval,probs=0.025),
    pval_97_5 = quantile(pval,probs=0.975)
  )


########################
#### Found Clusters ####
########################


res_file_loc<- paste("../summarized_results/pars_",par_no,"/clusters.csv",sep='')
res <-read.csv(res_file_loc)

summed_res <- res %>% group_by(hmax,rep) %>%
  summarise(hwcd = mean(hwcd),
            tp = sum(tp),
            tn = sum(tn),
            fp = sum(fp),
            fn = sum(fn),
            acc = (sum(tp)+sum(tn))/sum(tp,tn,fp,fn),
            prop_none = sum(n_found==0)/n(),
            prop_found = mean(n_found/n_true_clust),
            prop_exact_all = mean(n_exact/n_true_clust,na.rm=T),
            prop_exact_found = mean(n_exact/n_found,na.rm=T))



