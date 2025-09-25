###########################
#### Process dat frame ####
###########################

dat <- dat %>% mutate(
  ppv = tp/(tp+fp),
  tpr = tp/(tp+fn),
  f1 = 2 * ((tp/(tp+fp))*(tp/(tp+fn)) )/((tp/(tp+fp))+(tp/(tp+fn))),
  broad_coverage = n_broad /nrets,
  narrow_coverage= n_narrow/nrets,
  exact_coverage = n_exact /nrets,
  broad_mapping = broad_compat/est_rets,
  exact_mapping = exact_compat/est_rets
)

in_phy <- dat %>% group_by(phy) %>%
  summarise(weight= n(),
            across(where(is.numeric), ~mean(.x, na.rm = TRUE))
            ) %>%
  dplyr::select(-c(rep,hmax,
            tp,tn,fp,fn,
            est_rets,
            n_broad,n_narrow,n_exact,
            nrets,level,blobs,nt_blobs))
across_phy<- in_phy %>%
  summarise(across(where(is.numeric) & !all_of("weight"),
                   ~ weighted.mean(.x, w = weight, na.rm = TRUE)))
across_phy<-across_phy[,-1]

ppv_plot<-ggplot(in_phy,aes(x=ppv))+geom_histogram()+
  labs(title=file_name,x = "Positive Predictive Value", y = "Number of Networks")+
  theme_minimal()+xlim(-0.1,1.1)
ggsave(filename=paste(plot_dirs[1],file_name,'.png',sep=''),ppv_plot,device='png')

tpr_plot<-ggplot(in_phy,aes(x=tpr))+geom_histogram()+
  labs(title=file_name,x = "True Positive Rate", y = "Number of Networks")+
  theme_minimal()+xlim(-0.1,1.1)
ggsave(filename=paste(plot_dirs[2],file_name,'.png',sep=''),tpr_plot,device='png')

f1_plot<-ggplot(in_phy,aes(x=f1))+geom_histogram()+
  labs(title=file_name,x = "f1 statistic", y = "Number of Networks")+
  theme_minimal()+xlim(-0.1,1.1)
ggsave(filename=paste(plot_dirs[3],file_name,'.png',sep=''),f1_plot,device='png')

rm(ppv_plot,tpr_plot,f1_plot) #save memory, ride a cowboy

bro_cov_plot <- ggplot(in_phy,aes(x=broad_coverage))+geom_histogram()+
  labs(title='broad sense of found',x = "Average percent of true clusters with a compatible estimated cluster", y = "Number of Networks")+
  theme_minimal()+xlim(-0.1,1.1)
ggsave(filename=paste(plot_dirs[4],file_name,'.png',sep=''),bro_cov_plot,device='png')

ex_cov_plot <- ggplot(in_phy,aes(x=exact_coverage))+geom_histogram()+
  labs(title='exact sense of found',x = "Average percent of true clusters with a compatible estimated cluster", y = "Number of Networks")+
  theme_minimal()+xlim(-0.1,1.1)
ggsave(filename=paste(plot_dirs[5],file_name,'.png',sep=''),ex_cov_plot,device='png')

rm(bro_cov_plot,ex_cov_plot)

bro_map_plot <- ggplot(in_phy,aes(x=broad_mapping))+geom_histogram()+
  labs(title='broad sense',x = "Average percent of estimated clusters that map to a true cluster", y = "Number of Networks")+
  theme_minimal()+xlim(-0.1,1.1)
ggsave(filename=paste(plot_dirs[6],file_name,'.png',sep=''),bro_map_plot,device='png')

ex_map_plot <- ggplot(in_phy,aes(x=exact_mapping))+geom_histogram()+
  labs(title='exact sense of',x = "Average percent of estimated clusters that map to a true cluster", y = "Number of Networks")+
  theme_minimal()+xlim(-0.1,1.1)
ggsave(filename=paste(plot_dirs[7],file_name,'.png',sep=''),ex_map_plot,device='png')


rf_plot <- ggplot(in_phy,aes(x=rf_dist))+geom_histogram()+
  labs(x = "Average RF distance", y = "Number of Networks")+
  theme_minimal()
ggsave(filename=paste(plot_dirs[8],file_name,'.png',sep=''),rf_plot,device='png')

rm(bro_map_plot,ex_map_plot,rf_plot)

###############################
#### Process hyb_dat frame ####
###############################

hyb_props <- read.csv(paste('../data/pars_',par_no,'/hyb_dat.csv',sep=''))
hyb_dat <- merge(hyb_dat,hyb_props, by=c('phy','name'),all.x = T) %>%
  as_tibble() %>%
  group_by(phy,rep,blob_no) %>%
  mutate(n_in_blob = sum(exact)-exact) %>%
  ungroup()

in_phy_hyb <- hyb_dat %>% group_by(phy,name) %>%
  summarise(weight = n(),
            broad  = sum(broad),
            narrow = sum(narrow),
            exact  = sum(exact),
            broad_prop = sum(broad)/n(),
            exact_prop = sum(exact)/n()
            )
in_phy_hyb<- merge(in_phy_hyb,hyb_props, by=c('phy','name'),all.x = T)


nboots<- 1000
coef_dat<-data.frame(matrix(NA,nrow=nboots,ncol=4))
pval_dat<-data.frame(matrix(NA,nrow=nboots,ncol=4))

coef2_dat<-data.frame(matrix(NA,nrow=nboots,ncol=4))
pval2_dat<-data.frame(matrix(NA,nrow=nboots,ncol=4))

coef_dat3<-data.frame(matrix(NA,nrow=nboots,ncol=4))
pval_dat3<-data.frame(matrix(NA,nrow=nboots,ncol=4))

coef_dat4<-data.frame(matrix(NA,nrow=nboots,ncol=7))
pval_dat4<-data.frame(matrix(NA,nrow=nboots,ncol=7))

for(boot in 1:nboots){
  boot_dat <- hyb_dat %>% group_by(phy,name) %>% 
    slice_sample(n=1) %>% ungroup()
  my_form <- exact ~ dist50+external+blob_level+blob_size
  my_model <- summary(clogit(update(my_form, .~. +strata(phy)), data= boot_dat))
  pval_dat[boot,]<-my_model$coefficients[,5]
  coef_dat[boot,]<-my_model$coefficients[,1]
  
  my_form3 <- exact ~ dist50+cycled+blob_level+blob_size
  my_model3 <- summary(clogit(update(my_form3, .~. +strata(phy)), data= boot_dat))
  pval_dat3[boot,]<-my_model3$coefficients[,5]
  coef_dat3[boot,]<-my_model3$coefficients[,1]
  
  my_form4 <- exact ~ dist50+cy_adjacent+cycled+stacked+blob_level+blob_size+n_in_blob
  my_model4 <- summary(clogit(update(my_form4, .~. +strata(phy)), data= boot_dat))
  pval_dat4[boot,]<-my_model4$coefficients[,5]
  coef_dat4[boot,]<-my_model4$coefficients[,1]
  
  # start_vals <- c(0, 0, 0,0,0)
  # fit <- optim(par = start_vals,
  #              fn = full_negloglik,
  #              df = boot_dat,
  #              formula = my_form,
  #              method = "BFGS",
  #              hessian = TRUE,
  #              control = list(maxit = 1000))
  # 
  # my_model2<-summary.custom(fit)[-1,]
  # coef2_dat[boot,]<- my_model2[1,]
  # pval2_dat[boot,]<-my_model2[4,]
  
}
factor_names <-rownames(my_model$coefficients)
colnames(coef_dat)<-factor_names
colnames(pval_dat)<-factor_names
# colnames(coef2_dat)<-factor_names
# colnames(pval2_dat)<-factor_names

factor_names3 <-rownames(my_model3$coefficients)
colnames(coef_dat3)<-factor_names3
colnames(pval_dat3)<-factor_names3

factor_names4 <-rownames(my_model4$coefficients)
colnames(coef_dat4)<-factor_names4
colnames(pval_dat4)<-factor_names4

for(i in 1:ncol(coef_dat)){
  coef_name<-factor_names[i]
  pval_plot<-ggplot(data=pval_dat,aes(x=.data[[coef_name]]))+geom_histogram()+labs(x="p value")
  ggsave(filename=paste(plot_dirs[8+i],file_name,'.png',sep=''),pval_plot,device='png')
  
  coef_plot<- ggplot(data=coef_dat,aes(x=.data[[coef_name]]))+geom_histogram()+labs(x="Coefficient")
  ggsave(filename=paste(plot_dirs[12+i],file_name,'.png',sep=''),coef_plot,device='png')
  
  # pval2_plot<-ggplot(data=pval2_dat,aes(x=.data[[coef_name]]))+geom_histogram()+labs(x="p value")
  # ggsave(filename=paste(plot_dirs[8+i],file_name,'_custom.png',sep=''),pval2_plot,device='png')
  # 
  # coef2_plot<- ggplot(data=coef2_dat,aes(x=.data[[coef_name]]))+geom_histogram()+labs(x="Coefficient")
  # ggsave(filename=paste(plot_dirs[12+i],file_name,'_custom.png',sep=''),coef2_plot,device='png')
  # 
  
}
rm(pval_plot,coef_plot)

across_hyb <- in_phy_hyb %>%
  summarise(across(where(is.numeric) & !all_of("weight"),
                   ~ weighted.mean(.x, w = weight, na.rm = TRUE)))
across_hyb<-across_hyb[,-1]

###############################
#### Process tob_dat frame ####
###############################

in_phy_tob <- tob_dat %>% group_by(phy,rep) %>%
  summarise(est_nblobs = n(),
            n_compat = sum(true_blob != -1))
in_phy_tob <- in_phy_tob %>% group_by(phy) %>%
  summarise(weight=n(),
            perc_compat = mean(n_compat/est_nblobs))

across_tob <- in_phy_tob %>%
  summarise(across(where(is.numeric) & !all_of("weight"),
                   ~ weighted.mean(.x, w = weight, na.rm = TRUE)))

across_tob<-across_tob[,-1]

tob_plot<-ggplot(in_phy_tob,aes(x=perc_compat))+labs(x="Average percent of compatible blobs")+xlim(-0.1,1.1)+geom_histogram()
ggsave(filename=paste(plot_dirs[17],file_name,'.png',sep=''),plot=tob_plot,device='png')

rm(tob_plot)