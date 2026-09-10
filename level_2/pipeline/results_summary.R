library(tidyverse)
library(ggplot2)
library(survival)
library(xtable)
library(grafify)
library(patchwork)
library(ggridges)
library(vegan)



source("00.generate_seeds.R")
source("./functions.R")

setwd("../pipeline")



net_dat <- read.csv('../data/sim_net_props.csv')
net_dat <- net_dat %>%
  left_join(pars, by = c("par" = "setting_no"))

net_dat <- net_dat %>%
  # 1. Group the data by the specified variables
  group_by(nu, ntips, level1) %>%
  mutate(ID = cur_group_id()) %>%
  ungroup()

sim_sum_dat <- net_dat %>% group_by(ID) %>%
  summarise(nu = mean(nu),
            tips= mean(ntips),
            level1 = mean(level1),
            TC_ave = mean(TC),
            TB_ave = mean(TB),
            FU_ave =mean(FU),
            rets_ave =mean(nrets),
            rets_var =var(nrets),
            level_ave = mean(level),
            level_var = var(level))
sim_sum_dat$ID=NULL
xtable(sim_sum_dat)

rm(sim_sum_dat)
####
## plot nrets 
####
filt_dat <- net_dat%>% filter(level1==0)

my_labels <- c(
  "15" = "15 taxa", 
  "20" = "20 taxa", 
  "25" = "25 taxa"
)
save_dir <-'../figs/sim_sum/'
ggplot(filt_dat, aes(x = nrets, fill = as.factor(nu))) +
  geom_bar(position = "dodge",stat='count') +
  facet_wrap(~ ntips,labeller = labeller(ntips = my_labels)) + 
  theme_minimal() +
  labs(x='Number of reticulations', y='Frequency',fill='Hybridization rate')
ggsave(paste(save_dir,'nrets.png'),height=3,width=3)
ggplot(filt_dat, aes(x = level, fill = as.factor(nu))) +
  geom_bar(position = "dodge",stat='count') +
  facet_wrap(~ ntips,labeller = labeller(ntips = my_labels)) + 
  theme_minimal() +
  labs(x='Network level', y='Frequency',fill='Hybridization rate')
ggsave(paste(save_dir,'level.png'),height=3,width=3)
ggplot(filt_dat, aes(x = blobs, fill = as.factor(nu))) +
  geom_bar(position = "dodge",stat='count') +
  facet_wrap(~ ntips,labeller = labeller(ntips = my_labels)) + 
  theme_minimal() +
  labs(x='Number of blobs', y='Frequency',fill='Hybridization rate')
ggsave(paste(save_dir,'nblobs.png'),height=3,width=3)



rm(filt_dat)
####
#read in all branch length data
####


all_bls <-list()
for(par_no in 1:36){
  all_bls[[par_no]] <- read.csv(paste('../data/pars_',par_no,'/branch_lengths.csv',sep=''))
}
all_bls <- data.frame(
  setting_no = rep(seq_along(all_bls), sapply(all_bls, nrow)),
  branch_len = unlist(all_bls)
)

all_bls<-  all_bls  %>% left_join(pars,all_bls,by='setting_no') %>%
  mutate(ntips=as.factor(ntips),
         nu = as.factor(nu))
all_bls$ntips <- factor(all_bls$ntips, levels = c(25, 20, 15))

ggplot(data = all_bls, aes(x = branch_len, fill =ntips)) +
  geom_density(alpha = 0.7) +
  theme_minimal() +
  labs(y="Density of branch lengths",
       x = 'Branch length in coalescent units',
       fill = 'Number of taxa') +
  coord_cartesian(xlim = c(0, quantile(all_bls$branch_len, 0.99)))+
  scale_fill_grafify(palette = "fishy",reverse=TRUE)
  
ggsave(paste(save_dir,'branch_lengths.png',sep=''),height=3,width=4)
 
  
  #########################################
  #########################################
  #####                              ######
  #####       Summarize results      ######
  #####                              ######
  #########################################
  #########################################
  

all_dat<-list()
for(par_no in 1:18){
  print(par_no)
  res_file_loc <- paste("../summarized_results/pars_",par_no,"/",sep='')
  if(!file.exists(paste(res_file_loc,'clusters.csv',sep=''))){
    print(paste0('skipping ',par_no))
    next
  }

  cluster_res <- unique(read.csv(paste(res_file_loc,'clusters.csv',sep='')))
  tob_res     <-unique(read.csv(paste(res_file_loc,'tob.csv',sep='')))
  quar_res    <-unique(read.csv(paste(res_file_loc,'squirrel.csv',sep='')))


  
  tob_dat <- read.csv(paste(res_file_loc,'tob/compat.csv',sep=''))
  tob_dat <- tob_dat %>% group_by(phy,rep,hmax) %>%
    summarise(est_nblobs = n(),
              n_compat = sum(true_blob != -1), 
              n_blobs_found = length(unique(true_blob[true_blob != -1])))
  
  res<- merge(cluster_res,tob_res,by=c('phy','rep','hmax'),all=T)
  res<- merge(res,tob_dat,by=c('phy','rep','hmax'),all=T)
  res<- merge(res,quar_res,all=T)
  
  
  if(file.exists(paste(res_file_loc,'cf_dists.csv',sep=''))){
    dist_res    <-unique(read.csv(paste(res_file_loc,'cf_dists.csv',sep='')))
    dist_res_names <- c("phy","rep","hmax",
                        "CF_est_obs",
                        "CF_est_true",
                        "CF_obs_true",
                        "subnet_dist")
    colnames(dist_res)<-dist_res_names
    res<- merge(res,dist_res,by=c('phy','rep','hmax'),all=T)
    
  }

  par_dat <- net_dat %>% ##data from the true network
  filter(par==par_no) %>%
  dplyr::select(phy,nrets,level,nt_blobs,nu,level1,ntips,ngt) %>%
  rename(blobs = nt_blobs)
  
  all_dat[[par_no]] <- merge(res,par_dat,by=c('phy'),all.x=FALSE) #add true phy info
}
combined_df <- bind_rows(all_dat, .id = "setting_no")
combined_df <- combined_df %>%
  mutate(CF_dist = if_else(CF_dist < 0, NA, CF_dist),
         quar_found = if_else(quar_found < 0, NA, quar_found),
         quar_compat = if_else(quar_compat < 0, NA, quar_compat),
         CF_est_obs =  if_else(CF_est_obs < 0, NA, CF_est_obs),
         CF_est_true =  if_else(CF_est_true < 0, NA, CF_est_true),
         CF_obs_true =  if_else(CF_obs_true < 0, NA, CF_obs_true))
         
write.csv(combined_df,'../summarized_results/all_dat.csv')
rm(all_dat)



dat_filters = list(
  list(name='hmax1',filter = expr(filter(hmax==1))),
  list(name='hmax2',filter = expr(filter(hmax==2))),
  list(name='hmax3',filter = expr(filter(hmax==3))),
  list(name='hmax4',filter = expr(filter(hmax==4))),
  list(name='hmax5',filter = expr(filter(hmax==5))),
  list(name='blob',filter = expr(filter(pmin(blobs,max(hmax)) == hmax))),
  list(name='nrets',filter = expr(filter(pmin(nrets,max(hmax)) == hmax)))
)

diff_filters=list()
for(filt in dat_filters){
  diff_filters[[filt$name]] <- eval(expr(combined_df %>% group_by(setting_no,phy,rep)  %>% !!filt$filter %>% ungroup())) 
}

filter_df <- bind_rows(diff_filters, .id = "filter")
filter_df <- filter_df %>% mutate(
  ppv = tp/(tp+fp),
  tpr = tp/(tp+fn), #recall
  f1 = 2 * ((tp/(tp+fp))*(tp/(tp+fn)) )/((tp/(tp+fp))+(tp/(tp+fn))),
  broad_coverage = n_broad /nrets,
  narrow_coverage= n_narrow/nrets,
  exact_coverage = n_exact /nrets,
  broad_mapping = broad_compat/est_rets,
  exact_mapping = exact_compat/est_rets,
  perc_blob_compat = n_compat/est_nblobs,
  perc_blobs_found = n_blobs_found/blobs,
  found_sub = as.numeric(subnet_dist == 0)
) %>% dplyr::select( ##remove fields not needed
  -c(
     tp,tn,fp,fn,
     broad_compat,exact_compat)
)
write.csv(filter_df,'../summarized_results/filter_extended_dat.csv')

##take averages within reps and compute summary statistics
in_phy_filt <- filter_df %>% group_by(setting_no,phy,filter) %>%
  summarise(weight= n(),
            across(where(is.numeric), ~mean(.x, na.rm = TRUE))
  )  %>% dplyr::select( ##remove field now meaningless
    -c(rep)
  )

write.csv(in_phy_filt,'../summarized_results/rep_summed_filter_extended.csv')

rm(filter_df)


######################################
######################################
########                     #########
########   Excalidraw Plots  #########
########                     #########
######################################
######################################

in_phy_filt <- read.csv('../summarized_results/rep_summed_filter_extended.csv')


all_filts<-unique(in_phy_filt$filter)

for(filter_type in all_filts){
  print(paste("we are in filter",filter_type))
  # filter_type<-'nrets'
  save_dir <- paste('../figs/filter/',filter_type,'/',sep='')
  dir.create(save_dir,recursive=T)
  
  
  blob_phy_filt <- in_phy_filt %>%
    mutate(across(where(is.numeric), ~ if_else(.x < 0, NA, .x))) %>%
    filter(filter == filter_type) %>%
    mutate(
      # Convert grouping variables to factors so ggplot treats them as categories
      ntips_label = paste(ntips, "Tips"), # specific labels for the facets
      ngt = as.factor(ngt),
      level = if_else(level >= 5, "5≤", as.character(level)),
      level = factor(level, levels = c("1", "2", "3", "4", "5≤")),# bin
      # X-axis
    ) %>%
    ungroup()
  
  
  source('make_plots.R')
  
  
}





  
  
  ###############################
  ###############################
  ####                       ####
  ####  Plot across filters  ####
  ####                       ####
  ###############################
  ###############################
  ###############################








#########################
#### plot ppv,tpr,f1 ####
#########################

  

#filter by number of genes and level then plot
# ##Filter the data to just work with level2 data  that is 100 ngt

  
filter_type<-100 ##we summarize for a specific ngt 
in_phy_filt <- read.csv('../summarized_results/rep_summed_filter_extended.csv')



save_dir <- paste('../figs/filter_ngt/ngt_',filter_type,'/',sep='')
dir.create(save_dir,recursive=T)
  

##specify the blob and nrets as different from hmax
in_phy_filt$split_group <- ifelse(in_phy_filt$filter %in% c("blob", "nrets"), 
                                  "Controls","Hmax")
in_phy_filt$split_group <- factor(in_phy_filt$split_group,levels = c("Hmax", "Controls"))




line_data <- data.frame(split_group = "Hmax", x_pos = -Inf)
p1<-ggplot(in_phy_filt, aes(x = as.factor(filter), y = ppv, fill = as.factor(ntips))) +
  geom_boxplot(trim = F,
              position = position_dodge(0.7),width=1,
              draw_quantiles = 0.5) +
  geom_vline(data = line_data, aes(xintercept = x_pos), 
             linetype = "solid", 
             color = "grey50", 
             linewidth = 1) +
  facet_grid(~ split_group, scales = "free_x", space = "free_x") +
  theme_minimal()+
  theme(panel.spacing = unit(0.25, "lines"), # <--- Increases the gap size
        strip.text = element_blank()) +
  labs(y="Positive predictive value",
       x = "Hmax", 
       fill = "Taxa")+
  scale_fill_grafify(palette = "fishy",reverse=F)+
  coord_cartesian(clip = "off") +
  scale_x_discrete(labels = c("nrets" = "Reticulations","blob" = "Blob","hmax1" = "1","hmax2" = "2","hmax3" = "3","hmax4" = "4","hmax5" = "5"))


p2<-ggplot(in_phy_filt, aes(x = as.factor(filter), y = tpr, fill = as.factor(ntips))) +
  geom_boxplot(trim = F,
               position = position_dodge(0.7),width=1,
               draw_quantiles = 0.5) +
  geom_vline(data = line_data, aes(xintercept = x_pos), 
             linetype = "solid", 
             color = "grey50", 
             linewidth = 1) +
  facet_grid(~ split_group, scales = "free_x", space = "free_x") +
  theme_minimal()+
  theme(panel.spacing = unit(0.25, "lines"), # <--- Increases the gap size
        strip.text = element_blank()) +
  labs(y="True positive rate",
       x = "Hmax", 
       fill = "Taxa")+
  scale_fill_grafify(palette = "fishy",reverse=F)+
  coord_cartesian(clip = "off") +
  scale_x_discrete(labels = c("nrets" = "Reticulations","blob" = "Blob","hmax1" = "1","hmax2" = "2","hmax3" = "3","hmax4" = "4","hmax5" = "5"))


p3<-ggplot(in_phy_filt, aes(x = as.factor(filter), y = f1, fill = as.factor(ntips))) +
  geom_boxplot(trim = F,
               position = position_dodge(0.7),width=1,
               draw_quantiles = 0.5) +
  geom_vline(data = line_data, aes(xintercept = x_pos), 
             linetype = "solid", 
             color = "grey50", 
             linewidth = 1) +
  facet_grid(~ split_group, scales = "free_x", space = "free_x") +
  theme_minimal()+
  theme(panel.spacing = unit(0.25, "lines"), # <--- Increases the gap size
        strip.text = element_blank()) +
  labs(y="F1 statistic",
       x = "Hmax", 
       fill = "Taxa")+
  scale_fill_grafify(palette = "fishy",reverse=F)+
  coord_cartesian(clip = "off") +
  scale_x_discrete(labels = c("nrets" = "Reticulations","blob" = "Blob","hmax1" = "1","hmax2" = "2","hmax3" = "3","hmax4" = "4","hmax5" = "5"))


final_plot <- p1 + p2 + p3 + plot_layout(nrow = 1, guides = "collect") & theme(axis.text.x = element_text(angle = 45, hjust = 1))
final_plot
ggsave(paste(save_dir,'ppv_tpr.png',sep=''))

#########################
### plot hwcd #####
########################
clean_df <- in_phy_filt %>%
  ungroup() %>%  # <--- CRITICAL: Removes the list-column metadata
  mutate(
    ntips = as.factor(ntips),
    filter = factor(filter, levels = c('hmax1','hmax2','hmax3','hmax4','hmax5','nrets','blob')),
    # Create the split group logic again since this is a new df
    split_group = ifelse(filter %in% c("blob", "nrets"), "Controls", "Hmax"),
    split_group = factor(split_group, levels = c("Hmax", "Controls"))
  )

line_data <- data.frame(split_group = "Hmax", x_pos = -Inf)

ggplot(clean_df, aes(x = filter, y = hwcd, color = ntips)) +
  geom_line(data = subset(clean_df, split_group == "Hmax"),
            aes(group = interaction(phy, ntips,setting_no)), 
            alpha = 0.15) + 
  stat_summary(data = subset(clean_df, split_group == "Hmax"),
               aes(group = ntips),
               fun = mean, geom = "line", linewidth = 1.2) +
  geom_boxplot(data = subset(clean_df, split_group == "Controls"),
               aes(fill = ntips), 
               color = "grey30", alpha = 0.6, width = 0.6,
               position = position_dodge(0.8)) +
  geom_vline(data = line_data, aes(xintercept = x_pos), 
             linetype = "solid", color = "grey50") +
  
  facet_grid(~ split_group, scales = "free_x", space = "free_x") +
  
  scale_color_grafify(palette = "fishy") +
  scale_fill_grafify(palette = "fishy") +
  theme_minimal() +
  coord_cartesian(clip = "off") +
  theme(panel.spacing = unit(1, "lines"),
        strip.text = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  
  labs(y = "Hardwired Cluster Distance", x = "Hmax", color = "Taxa", fill = "Taxa") +
  scale_x_discrete(labels = c("nrets" = "Reticulations", "blob" = "Blob",
                              "hmax1" = "1", "hmax2" = "2", "hmax3" = "3",
                              "hmax4" = "4", "hmax5" = "5"))

ggsave(paste(save_dir,'hwcd_spaghetti.png',sep=''))



hwcd_data <- in_phy_filt %>% 
  filter(filter == 'blob' & nrets <= 7) 
tip_labels <- c("15" = "15 Tips", "20" = "20 Tips", "25" = "25 Tips")
ggplot(hwcd_data, aes(x = hwcd, y = as.factor(nrets), fill = as.factor(ntips))) +
  geom_density_ridges(alpha = 0.7, scale = 1.3, color = "white", size = 0.2) +
  facet_grid(ntips ~ ., 
             scales = "free_y", 
             space = "free_y", 
             labeller = as_labeller(tip_labels)) +
  

  scale_fill_grafify(palette = "fishy") +
  theme_minimal() +
  theme(
    strip.text.y = element_text(angle = 0, face = "bold", size = 10),
    strip.background = element_blank(),
    panel.spacing = unit(1, "lines"),
    panel.grid.major.y = element_line(color = "grey90", linetype = "dashed")
  ) +
  labs(title = "Hardwired Cluster Distance by Reticulations",
       y = "Number of Reticulations", 
       x = "Hardwired Cluster Distance (HWCD)",
       fill = "Taxa")

ggsave(paste(save_dir,'hwcd_ridge.png',sep=''))


######
##Coverage and mapping plots
######

plot_data <- in_phy_filt %>%
  pivot_longer(cols = c(broad_mapping, exact_mapping), 
               names_to = "mapping_type", 
               values_to = "mapping_value")
plot_data$split_group <- factor(ifelse(plot_data$filter %in% c("blob", "nrets"), 
                                       "Controls", "Hmax"), levels = c("Hmax", "Controls"))
ggplot(plot_data, aes(x = as.factor(filter), y = mapping_value, 
                      color = as.factor(ntips), 
                      group = interaction(as.factor(ntips), mapping_type))) + 
  stat_summary(data = subset(plot_data, split_group == "Hmax"),
               fun = mean, 
               geom = "line", 
               position=position_dodge(width = 0.4),
               aes(linetype = mapping_type), 
               linewidth = 1) +
  stat_summary(fun = mean, 
               geom = "point",
               position=position_dodge(width = 0.4),
               aes(shape = mapping_type),   
               size = 3) +
  stat_summary(fun.data = mean_se, 
               geom = "errorbar", 
               position=position_dodge(width = 0.4),
               width = 0.2) +
  geom_vline(data = line_data, aes(xintercept = x_pos), 
             linetype = "dashed", color = "grey50", linewidth = 1) +
  facet_grid(~ split_group, scales = "free_x", space = "free_x") +
  scale_color_grafify(palette = "fishy", reverse = FALSE) +
  theme_minimal() +
  coord_cartesian(clip = "off") +
  theme(panel.spacing = unit(1, "lines"), 
        strip.text = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Percent of hybridization mapped to the true network", x = "Hmax", 
       color = "Taxa", linetype = "Mapping", shape = "Mapping")+
  scale_linetype_discrete(labels = c("broad_mapping" = "Broad", 
                                     "exact_mapping" = "Exact")) +
  scale_shape_discrete(labels = c("broad_mapping" = "Broad", 
                                  "exact_mapping" = "Exact")) 
ggsave(paste(save_dir,'hyb_clad_map.png',sep=''))


plot_data <- in_phy_filt %>%
  pivot_longer(cols = c(narrow_coverage,broad_coverage, exact_coverage), 
               names_to = "mapping_type", 
               values_to = "mapping_value")
plot_data$split_group <- factor(ifelse(plot_data$filter %in% c("blob", "nrets"), 
                                       "Controls", "Hmax"), levels = c("Hmax", "Controls"))
ggplot(plot_data, aes(x = as.factor(filter), y = mapping_value, 
                      color = as.factor(ntips), 
                      group = interaction(as.factor(ntips), mapping_type))) + 
  stat_summary(data = subset(plot_data, split_group == "Hmax"),
               fun = mean, 
               geom = "line", 
               position=position_dodge(width = 0.4),
               aes(linetype = mapping_type), 
               linewidth = 1) +
  stat_summary(fun = mean, 
               geom = "point",
               position=position_dodge(width = 0.4),
               aes(shape = mapping_type),   
               size = 3) +
  #stat_summary(fun.data = mean_se, 
  #             geom = "errorbar", 
  #             position=position_dodge(width = 0.4),
  #             width = 0.2) +
  geom_vline(data = line_data, aes(xintercept = x_pos), 
             linetype = "dashed", color = "grey50", linewidth = 1) +
  facet_grid(~ split_group, scales = "free_x", space = "free_x") +
  scale_color_grafify(palette = "fishy", reverse = FALSE) +
  theme_minimal() +
  coord_cartesian(clip = "off") +
  theme(panel.spacing = unit(1, "lines"), 
        strip.text = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Percent of true hybridizations found in the estimated network", x = "Hmax", 
       color = "Taxa", linetype = "Coverage", shape = "Coverage")+
  scale_linetype_discrete(labels = c("broad_coverage" = "Broad", 
                                     "exact_coverage" = "Exact",
                                     "narrow_coverage" = "Narrow")) +
  scale_shape_discrete(labels = c("broad_coverage" = "Broad", 
                                  "exact_coverage" = "Exact",
                                  "narrow_coverage" = "Narrow")) 
ggsave(paste(save_dir,'hyb_clade_cov.png',sep=''))


##################################
####Plot Tob compat and found ####
##################################



in_phy_filt$split_group <- factor(ifelse(in_phy_filt$filter %in% c("blob", "nrets"), 
                                       "Controls", "Hmax"), levels = c("Hmax", "Controls"))
ggplot(in_phy_filt, aes(x = as.factor(filter), y = perc_blob_compat , 
                      fill = as.factor(ntips), 
                      )) +
  geom_boxplot(trim = F,
              position = position_dodge(0.7),width=1,
              draw_quantiles = 0.5) +
  geom_vline(data = line_data, aes(xintercept = x_pos), 
             linetype = "solid", 
             color = "grey50", 
             linewidth = 1) +
  facet_grid(~ split_group, scales = "free_x", space = "free_x") +
  theme_minimal()+
  theme(panel.spacing = unit(0.25, "lines"), # <--- Increases the gap size
        strip.text = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y="Estimated blobs that are compatible with a true blob",
       x = "Hmax", 
       fill = "Taxa")+
  scale_fill_grafify(palette = "fishy",reverse=F)+
  coord_cartesian(clip = "off") +
  scale_x_discrete(labels = c("nrets" = "Reticulations","blob" = "Blob","hmax1" = "1","hmax2" = "2","hmax3" = "3","hmax4" = "4","hmax5" = "5"))
ggsave(paste(save_dir,'blob_compat.png',sep=''))


ggplot(in_phy_filt, aes(x = as.factor(filter), y = perc_blobs_found, 
                        fill = as.factor(ntips), 
)) +
  geom_boxplot(trim = F,
               position = position_dodge(0.7),width=1,
               draw_quantiles = 0.5) +
  geom_vline(data = line_data, aes(xintercept = x_pos), 
             linetype = "solid", 
             color = "grey50", 
             linewidth = 1) +
  facet_grid(~ split_group, scales = "free_x", space = "free_x") +
  theme_minimal()+
  theme(panel.spacing = unit(0.25, "lines"), # <--- Increases the gap size
        strip.text = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y="True blob compatibility with estimated blobs",
       x = "Hmax", 
       fill = "Taxa")+
  scale_fill_grafify(palette = "fishy",reverse=F)+
  coord_cartesian(clip = "off") +
  scale_x_discrete(labels = c("nrets" = "Reticulations","blob" = "Blob","hmax1" = "1","hmax2" = "2","hmax3" = "3","hmax4" = "4","hmax5" = "5"))
ggsave(paste(save_dir,'blob_cov.png',sep=''))






#################
######################
# Plot ToB distances
#######################
#######################


plot_data <- in_phy_filt %>%
  pivot_longer(cols = c(rf_dist,split_dist), 
               names_to = "mapping_type", 
               values_to = "mapping_value")
plot_data$split_group <- factor(ifelse(plot_data$filter %in% c("blob", "nrets"), 
                                       "Controls", "Hmax"), levels = c("Hmax", "Controls"))
ggplot(plot_data, aes(x = as.factor(filter), y = mapping_value, 
                      color = as.factor(ntips), 
                      group = interaction(as.factor(ntips), mapping_type))) + 
  stat_summary(data = subset(plot_data, split_group == "Hmax"),
               fun = mean, 
               geom = "line", 
               position=position_dodge(width = 0.4),
               aes(linetype = mapping_type), 
               linewidth = 1) +
  stat_summary(fun = mean, 
               geom = "point",
               position=position_dodge(width = 0.4),
               aes(shape = mapping_type),   
               size = 3) +
  #stat_summary(fun.data = mean_se, 
  #             geom = "errorbar", 
  #             position=position_dodge(width = 0.4),
  #             width = 0.2) +
  geom_vline(data = line_data, aes(xintercept = x_pos), 
             linetype = "dashed", color = "grey50", linewidth = 1) +
  facet_grid(~ split_group, scales = "free_x", space = "free_x") +
  scale_color_grafify(palette = "fishy", reverse = FALSE) +
  theme_minimal() +
  coord_cartesian(clip = "off") +
  theme(panel.spacing = unit(1, "lines"), 
        strip.text = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Distance between tree of blobs", x = "Hmax", 
       color = "Taxa", linetype = "Distance measure", shape = "Distance measure")+
  scale_linetype_discrete(labels = c("rf_dist" = "Robinson-Foulds", 
                                     "split_dist" = "Matching Split"))+
  scale_shape_discrete(labels = c("rf_dist" = "Robinson-Foulds", 
                                  "split_dist" = "Matching Split")) 
ggsave(paste(save_dir,'tob_dists.png',sep=''))


ridge_data <- in_phy_filt %>% 
  filter(filter == 'blob' & nrets <= 7) %>%
  pivot_longer(cols = c(split_dist, rf_dist), 
               names_to = "metric_name", 
               values_to = "dist_value")

tip_labels <- c("15" = "15 Tips", "20" = "20 Tips", "25" = "25 Tips")

ggplot(ridge_data, aes(x = dist_value, y = as.factor(nrets), fill = metric_name)) +
  geom_density_ridges(alpha = 0.6, scale = 1.2, color = "white", size = 0.2) +
  facet_grid(ntips ~ ., 
             scales = "free_y", 
             space = "free_y", 
             labeller = as_labeller(tip_labels)) +
  scale_fill_manual(
    values = c("split_dist" = "#E69F00", "rf_dist" = "#56B4E9"),
    labels = c("split_dist" = "Matching split", "rf_dist" = "Robinson-Foulds") 
  ) +
  theme_minimal() + 
  theme(
    strip.text.y = element_text(angle = 0, face = "bold", size = 10),
    strip.background = element_blank(),
    panel.spacing = unit(1, "lines")
  ) +
  labs( 
       fill = "Metric",
       y = "Number of Reticulations", 
       x = "Distance") +
  xlim(0, 10)

ggsave(paste(save_dir,'tob_dists_ridge.png',sep=''))





####


hyb_res <- read.csv(paste("../summarized_results/pars_",par_no,"/found_clades.csv",sep=''))
hyb_res <- merge(hyb_res,par_dat,by='phy',all.x=TRUE) #add true phy info

tob_res <- read.csv(paste("../summarized_results/pars_",par_no,"/tob/compat.csv",sep=''))
tob_res <- merge(tob_res,par_dat,by='phy',all.x=TRUE) #add true phy info


#############
#############
  

dat_filters = list(
  list(name='hmax1',filter = expr(filter(hmax==1))),
  list(name='hmax2',filter = expr(filter(hmax==2))),
  list(name='hmax3',filter = expr(filter(hmax==3))),
  list(name='hmax4',filter = expr(filter(hmax==4))),
  list(name='hmax5',filter = expr(filter(hmax==5))),
  list(name='blob',filter = expr(filter(min(blobs,5) == hmax))),
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



