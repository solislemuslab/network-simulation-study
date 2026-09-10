library(tidyverse)
library(ggplot2)
library(survival)
library(xtable)
library(grafify)
library(patchwork)
library(ggridges)
library(vegan)
library(broom) # Added broom to easily extract model coefficients cleanly
library(foreach)
library(doParallel)
library(rpart)
library(caret)
library(ipred)

source("00.generate_seeds.R")
source("./functions.R")

setwd("../pipeline")


net_dat <- read.csv('../data/sim_net_props.csv')
net_dat <- net_dat %>%
  left_join(pars, by = c("par" = "setting_no"))


all_dat<-list()
for(par_no in 1:36){
  
  hyb_dat_loc <- paste("../summarized_results/pars_",par_no,'/found_clades.csv',sep='')
  hyb_props_loc <- paste('../data/pars_',par_no,'/hyb_dat.csv',sep='')
  if(!file.exists(hyb_props_loc)){
    print(paste("no sim hyb dat for par no",par_no))
    next
  }
  if(!file.exists(hyb_dat_loc)){
    next
  }
  
  
  hyb_dat <- read.csv(hyb_dat_loc)
  hyb_props <- read.csv(hyb_props_loc)
  
  
  hyb_dat <- merge(hyb_dat,hyb_props, by=c('phy','name'),all.x = T) %>%
    as_tibble() %>%
    group_by(phy,rep,blob_no) %>%
    mutate(n_in_blob = sum(exact)-exact) %>%
    ungroup()
  
  all_dat[[par_no]] <- hyb_dat
}

combined_df <- all_dat %>%
  keep(~ nrow(.x) > 0) %>%
  bind_rows(.id = "setting_no") %>%
  # Create a globally unique network ID so strata doesn't mix different settings
  mutate(unique_phy = paste(setting_no, phy, sep = "_"))
rm(all_dat)


n_boots <- 100 

# 1. Setup the parallel backend
# Use detectCores() - 1 to leave one core free for your operating system
num_cores <- parallel::detectCores() - 2
# outfile = "" ensures that your print() statements will still show up in the console
cl <- makeCluster(num_cores, outfile = "")
registerDoParallel(cl)

###Run bootstraps in parallel
all_boot_results <- foreach(i = 1:n_boots, .combine = dplyr::bind_rows, .packages = c("dplyr", "survival", "broom", "caret", "rpart", "ipred")) %dopar% {
  
  # Print every 10 iterations
  if (i %% 10 == 0) {
    cat(sprintf("Heartbeat: Starting bootstrap iteration %d of %d...\n", i, n_boots))
  }
  
  # Ensure reproducible random sampling across parallel workers
  set.seed(42 + i) 
  

  my_form_clogit <- exact ~ dist50+external+blob_level+blob_size+cycle_size+cycled+cy_adjacent+n_in_blob+hmax
  my_form_bag <- exact_factor ~ dist50+external+blob_level+blob_size+cycle_size+cycled+cy_adjacent+n_in_blob+hmax
  
  #Slice the data for 1 rep
  boot_dat <- combined_df %>% 
    distinct(setting_no, phy, hmax, rep) %>% 
    group_by(setting_no, phy, hmax) %>% 
    slice_sample(n=1) %>% 
    ungroup() %>%
    inner_join(combined_df, by = c("setting_no", "phy", "hmax", "rep")) %>%
    # Create factor target variable for Caret classification
    mutate(exact_factor = factor(exact, levels = c(0, 1), labels = c("No", "Yes")))
  
  ###########################################
  #######Conditional Logistic Regression#####
  ###########################################
  
  model <- clogit(update(my_form_clogit, .~. +strata(unique_phy)), data= boot_dat)
  clogit_res <- tidy(model) %>%
    mutate(boot_iter = i, analysis = "clogit")
  
  ########################################
  #####Bagged Decision Trees (Caret)######
  ########################################
  
  # Partition by the outcome factor to maintain class balance in Train/Test sets
  train_index <- caret::createDataPartition(boot_dat$exact_factor, p = 0.8, list = FALSE)
  train_data <- boot_dat[train_index, ]
  test_data <- boot_dat[-train_index, ]
  
  ctrl <- caret::trainControl(method = "cv", number = 10) 
  
  bagged_cv <- caret::train(
    my_form_bag,
    data = train_data,
    method = "treebag",
    trControl = ctrl,
    importance = TRUE
  )
  
  # Extract Variable Importance
  imp <- caret::varImp(bagged_cv)$importance
  bag_imp_res <- data.frame(term = rownames(imp), estimate = imp$Overall) %>%
    mutate(boot_iter = i, analysis = "treebag_importance")
  
  # Extract Test Accuracy
  preds <- predict(bagged_cv, newdata = test_data)
  acc <- mean(preds == test_data$exact_factor, na.rm = TRUE)
  bag_acc_res <- data.frame(term = "Test_Accuracy", estimate = acc) %>%
    mutate(boot_iter = i, analysis = "treebag_accuracy")
  
  # Return combined results for this bootstrap iteration
  bind_rows(clogit_res, bag_imp_res, bag_acc_res)
}

stopCluster(cl)

#############################################
#### SUMMARIZE RESULTS ACROSS BOOTSTRAPS ####
#############################################
write.csv(all_boot_results,"../summarized_results/boot_regression.csv")
all_boot_results <-read.csv("../summarized_results/boot_regression.csv")


#Summarize CLOGIT Coefficients
boot_summary_clogit <- all_boot_results %>%
  filter(analysis == "clogit") %>%
  group_by(term) %>%
  summarise(
    mean_estimate = mean(estimate),
    boot_se = sd(estimate),
    ci_lower = quantile(estimate, 0.025),
    ci_upper = quantile(estimate, 0.975),
    z_stat = mean(estimate) / sd(estimate),
    p_value = 2 * pnorm(-abs(mean(estimate) / sd(estimate)))
  ) %>%
#  mutate(
#    Odds_Ratio = exp(mean_estimate),
#    OR_ci_lower = exp(ci_lower),
#    OR_ci_upper = exp(ci_upper)
#  ) %>%
  arrange(desc(p_value))

print("--- CLOGIT SUMMARY ---")
print(boot_summary_clogit)

xtable(boot_summary_clogit)

#Summarize TREEBAG Variable Importance
boot_summary_treebag_imp <- all_boot_results %>%
  filter(analysis == "treebag_importance") %>%
  group_by(term) %>%
  summarise(
    mean_importance = mean(estimate),
    sd_importance = sd(estimate),
    ci_lower = quantile(estimate, 0.025),
    ci_upper = quantile(estimate, 0.975)
  ) %>%
  arrange(desc(mean_importance))

print("--- BAGGED TREE VARIABLE IMPORTANCE SUMMARY ---")
print(boot_summary_treebag_imp)

xtable(boot_summary_treebag_imp)

#Summarize TREEBAG Accuracy
boot_summary_treebag_acc <- all_boot_results %>%
  filter(analysis == "treebag_accuracy") %>%
  summarise(
    mean_accuracy = mean(estimate),
    sd_accuracy = sd(estimate),
    ci_lower = quantile(estimate, 0.025),
    ci_upper = quantile(estimate, 0.975)
  )

print(boot_summary_treebag_acc)
