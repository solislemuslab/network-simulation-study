
#######################
##### CF distance #####
#######################

ggplot(blob_phy_filt, aes(x = level, y = CF_dist, fill = ngt)) +
  geom_violin(position = position_dodge(width = 0.8), 
              alpha = 0.6,          # Slightly more transparent to see gridlines
              trim = FALSE, 
              scale = "width",
              color = "grey50",     
              size = 0.3) +        
  geom_boxplot(position = position_dodge(width = 0.8), 
               width = 0.15,        
               color = "grey20",    
               fill = "white",      
               alpha = 0.9,         
               outlier.shape = NA) +
  facet_wrap(~ ntips_label) +
  scale_fill_grafify(palette = "fishy") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(), # Cleaner background
    strip.text = element_text(size = 12, face = "bold"),
    axis.text = element_text(color = "grey30"),
    plot.title = element_text(face = "bold", size = 14)
  ) +
  labs(
    x = "Network Level",
    y = expression(group("|", CF[obs] - CF[exp], "|")),
    fill = "Number of Gene Trees"
  )
ggsave(paste(save_dir,'CF_dist.png',sep=''))

#############
### GOF #####
#############

ggplot(blob_phy_filt, aes(x = level, y = gof_pval, fill = ngt)) +
  geom_violin(position = position_dodge(width = 0.8), 
              alpha = 0.6,          # Slightly more transparent to see gridlines
              trim = T, 
              scale = "width",
              color = "grey50",     
              size = 0.3) +        
  geom_boxplot(position = position_dodge(width = 0.8), 
               width = 0.15,        
               color = "grey20",    
               fill = "white",      
               alpha = 0.9,         
               outlier.shape = NA) +
  facet_wrap(~ ntips_label) +
  scale_fill_grafify(palette = "fishy") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(), # Cleaner background
    strip.text = element_text(size = 12, face = "bold"),
    axis.text = element_text(color = "grey30"),
    plot.title = element_text(face = "bold", size = 14)
  ) +
  labs(
    x = "Network Level",
    y = "GoF p-value",
    fill = "Number of Gene Trees"
  )
ggsave(paste(save_dir,'GoF_violin.png',sep=''))


tip_labels <- c("15" = "15 Tips", "20" = "20 Tips", "25" = "25 Tips")
ggplot(blob_phy_filt, aes(x = gof_pval, y = level, fill = as.factor(ntips))) +
  geom_density_ridges(alpha = 0.7, scale = 1.3, color = "white", size = 0.2) +
  facet_grid(ntips ~ ., 
             scales = "free_y", 
             space = "free_y", 
             labeller = as_labeller(tip_labels)) +
  scale_fill_grafify(palette = "fishy") + # Or use scale_fill_viridis_d()
  theme_minimal() +
  
  theme(
    strip.text.y = element_text(angle = 0, face = "bold", size = 10),
    strip.background = element_blank(),
    panel.spacing = unit(1, "lines"),
    panel.grid.major.y = element_line(color = "grey90", linetype = "dashed")
  ) +
  
  labs(title = "Quartet Goodness-of-Fit",
       y = "Number of Reticulations", 
       x = "p-value",
       fill = "Taxa")+
  xlim(c(0,1))
ggsave(paste(save_dir,'GoF_ridge.png',sep=''))

####################
#### TPR vs PPV ####
######################


cols <- c("Recall" = "#1A5F7A", "Precision" = "#D86C4E")
ggplot(blob_phy_filt, aes(x = level, group = ngt, linetype = ngt)) +
  
  # --- LEFT AXIS DATA (RECALL / TPR) ---
  stat_summary(aes(y = tpr, color = "Recall"), 
               fun = mean, geom = "line", size = 1.2) +
  stat_summary(aes(y = tpr, color = "Recall"), 
               fun = mean, geom = "point", size = 3) +
  # --- RIGHT AXIS DATA (PRECISION / PPV) ---
  stat_summary(aes(y = ppv, color = "Precision"), 
               fun = mean, geom = "line", size = 1.2) +
  stat_summary(aes(y = ppv, color = "Precision"), 
               fun = mean, geom = "point", size = 3) +
  # --- DUAL AXIS SETUP ---
  scale_y_continuous(
    name = "Recall (TPR)",
    sec.axis = sec_axis(~ ., name = "Precision (PPV)")
  ) +
  facet_wrap(~ ntips_label) +
  scale_color_manual(values = cols) + ##manual color wrapping
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 12, face = "bold"),
    axis.title.y.left = element_text(color = cols["Recall"], face = "bold"),
    axis.title.y.right = element_text(color = cols["Precision"], face = "bold")
  ) +
  labs(
    x = "Network Level",
    color = "Metric",
    linetype = "Gene Trees")
ggsave(paste(save_dir,'recall_precision_hm.png',sep=''))


#############################
##### HWCD Canon and FU dists #####
#############################


cols <- c("Canonical HWCD" = "#1A5F7A", "FU-Stable HWCD" = "#D86C4E")
ggplot(blob_phy_filt, aes(x = level, group = ngt, linetype = ngt)) +
  
  # --- LEFT AXIS DATA 
  stat_summary(aes(y = canon_hwcd, color = "Canonical HWCD"), 
               fun = mean, geom = "line", size = 1.2) +
  stat_summary(aes(y = canon_hwcd, color = "Canonical HWCD"), 
               fun = mean, geom = "point", size = 3) +
  # --- RIGHT AXIS DATA 
  stat_summary(aes(y = fu_hwcd, color = "FU-Stable HWCD"), 
               fun = mean, geom = "line", size = 1.2) +
  stat_summary(aes(y = fu_hwcd, color = "FU-Stable HWCD"), 
               fun = mean, geom = "point", size = 3) +
  # --- DUAL AXIS SETUP ---
  scale_y_continuous(
    name = "Canonical HWCD",
    sec.axis = sec_axis(~ ., name = "FU-Stable HWCD")
  ) +
  facet_wrap(~ ntips_label) +
  scale_color_manual(values = cols) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 12, face = "bold"),
    axis.title.y.left = element_text(color = cols["Canonical HWCD"], face = "bold"),
    axis.title.y.right = element_text(color = cols["FU-Stable HWCD"], face = "bold")
  ) +
  labs(
    x = "Network Level",
    color = "Metric",
    linetype = "Gene Trees")
ggsave(paste(save_dir,'reduced_dists.png',sep=''))


ggplot(blob_phy_filt, aes(x = level, y = hwcd, fill = ngt)) +
  geom_violin(position = position_dodge(width = 0.8), 
              alpha = 0.6,          # Slightly more transparent to see gridlines
              trim = FALSE, 
              scale = "width",
              color = "grey50",     
              size = 0.3) +        
  geom_boxplot(position = position_dodge(width = 0.8), 
               width = 0.15,        
               color = "grey20",    
               fill = "white",      
               alpha = 0.9,         
               outlier.shape = NA) +
  facet_wrap(~ ntips_label) +
  scale_fill_grafify(palette = "fishy") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(), # Cleaner background
    strip.text = element_text(size = 12, face = "bold"),
    axis.text = element_text(color = "grey30"),
    plot.title = element_text(face = "bold", size = 14)
  ) +
  labs(
    x = "Network Level",
    y = 'Hardwired Cluster Distance',
    fill = "Number of Gene Trees"
  )
ggsave(paste(save_dir,'HWCD.png',sep=''))


scale_factor <- 0.1
cols <- c("total" = "#E69F00", "average" = "#56B4E9")

ggplot(blob_phy_filt, aes(x = level, group = ngt, linetype = ngt)) +
  
  stat_summary(aes(y = total_disp_rf, color = "total"), 
               fun = mean, geom = "line", size = 1.2) +
  stat_summary(aes(y = total_disp_rf, color = "total"), 
               fun = mean, geom = "point", size = 3) +
  
  stat_summary(aes(y = ave_disp_rf / scale_factor, color = "average"), 
               fun = mean, geom = "line", size = 1.2) +
  stat_summary(aes(y = ave_disp_rf / scale_factor, color = "average"), 
               fun = mean, geom = "point", size = 3) +
  
  # --- DUAL AXIS SETUP ---
  scale_y_continuous(
    name = "Total paired RF distance",
    # We MULTIPLY by scale_factor here to make the labels show the comparable numbers
    sec.axis = sec_axis(~ . * scale_factor, name = "Average paired RF")
  ) +
  facet_wrap(~ ntips_label) +
  scale_color_manual(values = cols) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 12, face = "bold"),
    # Match Axis Colors for Clarity
    axis.title.y.left = element_text(color = cols["total"], face = "bold"),
    axis.title.y.right = element_text(color = cols["average"], face = "bold")
  ) +
  labs(x = "Network Level",
       color = "Metric",
       linetype = "Gene Trees")
ggsave(paste(save_dir,'displayed_hwcd.png',sep=''))





###################
## TOB Distances ##
###################

scale_factor <- 1

cols <- c("RF" = "#E69F00", "SplitInfo" = "#56B4E9")

ggplot(blob_phy_filt, aes(x = level, group = ngt, linetype = ngt)) +
  
  # --- LEFT AXIS: RF DISTANCE (Primary) ---
  stat_summary(aes(y = rf_dist, color = "RF"), 
               fun = mean, geom = "line", size = 1.2) +
  stat_summary(aes(y = rf_dist, color = "RF"), 
               fun = mean, geom = "point", size = 3) +
  # --- RIGHT AXIS: SPLIT INFO (Secondary - Scaled Down) ---
  # We DIVIDE by scale_factor here to bring the line down visually
  stat_summary(aes(y = split_dist / scale_factor, color = "SplitInfo"), 
               fun = mean, geom = "line", size = 1.2) +
  stat_summary(aes(y = split_dist / scale_factor, color = "SplitInfo"), 
               fun = mean, geom = "point", size = 3) +
  # --- DUAL AXIS SETUP ---
  scale_y_continuous(
    name = "Robinson-Foulds Distance",
    # We MULTIPLY by scale_factor here to make the labels show the comparable numbers
    sec.axis = sec_axis(~ . * scale_factor, name = "Split Information Distance")
  ) +
  facet_wrap(~ ntips_label) +
  scale_color_manual(values = cols) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 12, face = "bold"),
    axis.title.y.left = element_text(color = cols["RF"], face = "bold"),
    axis.title.y.right = element_text(color = cols["SplitInfo"], face = "bold")
  ) +
  labs(title = "Topological Accuracy: RF vs Split Information Distance",
       x = "Network Level",
       color = "Metric",
       linetype = "Gene Trees")
ggsave(paste(save_dir,'ToB_dists.png',sep=''))

############################################
####### Quarnet consistency scores #########
############################################


cols <- c("Found Quarnets" = "#1A5F7A", "Compatible Quarnets" = "#D86C4E")
ggplot(blob_phy_filt, aes(x = level, group = ngt, linetype = ngt)) +
  
  # --- LEFT AXIS DATA 
  stat_summary(aes(y = quar_found, color = "Found Quarnets"), 
               fun = mean, geom = "line", size = 1.2) +
  stat_summary(aes(y = quar_found, color = "Found Quarnets"), 
               fun = mean, geom = "point", size = 3) +
  # --- RIGHT AXIS DATA 
  stat_summary(aes(y = quar_compat, color = "Compatible Quarnets"), 
               fun = mean, geom = "line", size = 1.2) +
  stat_summary(aes(y = quar_compat, color = "Compatible Quarnets"), 
               fun = mean, geom = "point", size = 3) +
  # --- DUAL AXIS SETUP ---
  scale_y_continuous(
    name = "Found Quarnets",
    sec.axis = sec_axis(~ ., name = "Compatible Quarnets")
  ) +
  facet_wrap(~ ntips_label) +
  scale_color_manual(values = cols) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 12, face = "bold"),
    axis.title.y.left = element_text(color = cols["Found Quarnets"], face = "bold"),
    axis.title.y.right = element_text(color = cols["Compatible Quarnets"], face = "bold")
  ) +
  labs( title = 'Quarnet Compatability Score',
        x = "Network Level",
        color = "Metric",
        linetype = "Gene Trees")
ggsave(paste(save_dir,'quarnet_score.png',sep=''))



##############################
##############################
##### Hybrid clade heatmap ###
##############################
##############################

##### exact coverage ######

plot_dat <- blob_phy_filt %>%  group_by(ntips, level, ngt) %>%
  summarise(mean_pct = mean(exact_coverage, na.rm = TRUE) * 100, .groups = "drop") %>%
  mutate(
    ntips_label = paste(ntips, "Tips"),
    ngt = as.factor(ngt),
    level = as.factor(level)
  )


ggplot(plot_dat, aes(x = ngt, y = level, fill = mean_pct)) +
  geom_tile(color = "white", lwd = 0.5) + # White borders make it look cleaner
  geom_text(aes(label = round(mean_pct, 1)), 
            color = "white", size = 3, fontface = "bold") +

  facet_wrap(~ ntips_label) +
  scale_fill_gradientn(colors = c("#440154", "#21908C", "#FDE725"), 
                       name = "% Found") +
  
  theme_minimal() +
  coord_fixed() + # Keeps the tiles square
  theme(
    panel.grid = element_blank(), # Heatmaps don't need gridlines
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "right"
  ) +
  labs(title = "Hybrid Cluster Exact Coverage",
       subtitle = "Percentage of true clusters exactly found in estimated network",
       x = "Number of gene trees",
       y = "Network level")
ggsave(paste(save_dir,'hyb_clade_ex_cov_hm.png',sep=''))


######exact mapping #######

plot_dat <- blob_phy_filt %>%  group_by(ntips, level, ngt) %>%
  summarise(mean_pct = mean(exact_mapping, na.rm = TRUE) * 100, .groups = "drop") %>%
  mutate(
    ntips_label = paste(ntips, "Tips"),
    ngt = as.factor(ngt),
    level = as.factor(level)
  )
ggplot(plot_dat, aes(x = ngt, y = level, fill = mean_pct)) +
  geom_tile(color = "white", lwd = 0.5) +
  geom_text(aes(label = round(mean_pct, 1)), 
            color = "white", size = 3, fontface = "bold") +
  facet_wrap(~ ntips_label) +
  scale_fill_gradientn(colors = c("#440154", "#21908C", "#FDE725"), 
                       name = "% compatible") +
  theme_minimal() +
  coord_fixed() + # Keeps the tiles square
  theme(
    panel.grid = element_blank(), # Heatmaps don't need gridlines
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "right"
  ) +
  labs(title = "Hybrid Cluster Exact Mapping",
       subtitle = "Percentage of estimated clusters that exactly match the true network",
       x = "Number of gene trees",
       y = "Network level")
ggsave(paste(save_dir,'hyb_clade_ex_mapping_hm.png',sep=''))


###### broad mapping

plot_dat <- blob_phy_filt %>%  group_by(ntips, level, ngt) %>%
  summarise(mean_pct = mean(broad_mapping, na.rm = TRUE) * 100, .groups = "drop") %>%
  mutate(
    ntips_label = paste(ntips, "Tips"),
    ngt = as.factor(ngt),
    level = as.factor(level)
  )

ggplot(plot_dat, aes(x = ngt, y = level, fill = mean_pct)) +
  geom_tile(color = "white", lwd = 0.5) + # White borders make it look cleaner
  geom_text(aes(label = round(mean_pct, 1)), 
            color = "white", size = 3, fontface = "bold") +
  facet_wrap(~ ntips_label) +
  scale_fill_gradientn(colors = c("#440154", "#21908C", "#FDE725"), 
                       name = "% compatible") +
  theme_minimal() +
  coord_fixed() + # Keeps the tiles square
  theme(
    panel.grid = element_blank(), # Heatmaps don't need gridlines
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "right"
  ) +
  labs(title = "Hybrid Cluster Broad Mapping",
       subtitle = "Percentage of estimated clusters that are compatible the true network",
       x = "Number of gene trees",
       y = "Network level")
ggsave(paste(save_dir,'hyb_clade_broad_map_hm.png',sep=''))


###broad coverage


plot_dat <- blob_phy_filt %>%  group_by(ntips, level, ngt) %>%
  summarise(mean_pct = mean(broad_coverage, na.rm = TRUE) * 100, .groups = "drop") %>%
  mutate(
    ntips_label = paste(ntips, "Tips"),
    ngt = as.factor(ngt),
    level = as.factor(level)
  )
ggplot(plot_dat, aes(x = ngt, y = level, fill = mean_pct)) +
  geom_tile(color = "white", lwd = 0.5) +
  geom_text(aes(label = round(mean_pct, 1)), 
            color = "white", size = 3, fontface = "bold") +
  facet_wrap(~ ntips_label) +
  scale_fill_gradientn(colors = c("#440154", "#21908C", "#FDE725"), 
                       name = "% Found") +
  theme_minimal() +
  coord_fixed() +
  theme(
    panel.grid = element_blank(), # Heatmaps don't need gridlines
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "right"
  ) +
  labs(title = "Hybrid Cluster Broad Coverage",
       subtitle = "Percentage of true clusters that are compatibile with the estimated network",
       x = "Number of gene trees",
       y = "Network level")
ggsave(paste(save_dir,'hyb_clade_broad_cov_hm.png',sep=''))


##############################
##############################
#####   ToB heatmap        ###
##############################
##############################

##### compatible ######

##make percentages 
plot_dat <- blob_phy_filt %>%  group_by(ntips, level, ngt) %>%
  summarise(mean_pct = mean(perc_blob_compat, na.rm = TRUE) * 100, .groups = "drop") %>%
  mutate(
    ntips_label = paste(ntips, "Tips"),
    ngt = as.factor(ngt),
  )



ggplot(plot_dat, aes(x = ngt, y = level, fill = mean_pct)) +
  geom_tile(color = "white", lwd = 0.5) + # White borders make it look cleaner
  geom_text(aes(label = round(mean_pct, 1)), 
            color = "white", size = 3, fontface = "bold") +
  facet_wrap(~ ntips_label) +
  scale_fill_gradientn(colors = c("#440154", "#21908C", "#FDE725"), 
                       name = "% Compatible") +
  theme_minimal() +
  coord_fixed() + # Keeps the tiles square
  theme(
    panel.grid = element_blank(), # Heatmaps don't need gridlines
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "right"
  ) +
  labs(title = "Blob Compatability",
       subtitle = "Percentage of blobs in the estimated network compatible with the true network",
       x = "Number of gene trees",
       y = "Network level")
ggsave(paste(save_dir,'tob_compat_hm.png',sep=''))


##### found ######

plot_dat <- blob_phy_filt %>%  group_by(ntips, level, ngt) %>%
  summarise(mean_pct = mean(perc_blobs_found, na.rm = TRUE) * 100, .groups = "drop") %>%
  mutate(
    ntips_label = paste(ntips, "Tips"),
    ngt = as.factor(ngt),
  )


ggplot(plot_dat, aes(x = ngt, y = level, fill = mean_pct)) +
  geom_tile(color = "white", lwd = 0.5) + # White borders make it look cleaner
  geom_text(aes(label = round(mean_pct, 1)), 
            color = "white", size = 3, fontface = "bold") +
  facet_wrap(~ ntips_label) +
  scale_fill_gradientn(colors = c("#440154", "#21908C", "#FDE725"), 
                       name = "% Found") +
  theme_minimal() +
  coord_fixed() + # Keeps the tiles square
  theme(
    panel.grid = element_blank(), # Heatmaps don't need gridlines
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "right"
  ) +
  labs(title = "Found Blobs",
       subtitle = "Percentage true blobs that were found in the estimated network",
       x = "Number of gene trees",
       y = "Network level")
ggsave(paste(save_dir,'tob_found_hm.png',sep=''))


