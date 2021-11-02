library(TreeSim)

set.seed(1453)

# brlengths around the root
brlengths <- c(0.050000000, 0.100000000, 0.500000000, 1.000000000, 2.000000000)

# simulate a tree
tree <- sim.bd.taxa(n = 30, numbsim = 1, lambda = 1, mu = 0, frac = 1, complete = TRUE, stochsampling = FALSE)

# drop the tips that will create a long branch for the root
drop_group <- c("t30", "t7", "t13", "t15", "t28", "t14", "t11", "t22", "t16", "t26", "t17", "t24", "t4", "t3", "t19", "t21", "t29", "t25")
cropped_tree <- drop.tip(tree[[1]], tip = drop_group)

# write cropped tree to backup newick
write.tree(phy = cropped_tree, file = "cropped_backup_tree.newick")

# write to str object
cropped_tree_str <- write.tree(phy = cropped_tree)

# declare objects related to brlengths

short <- 1.077761864
long <- 2.029328697
diff_brlength <- long - short

# use placeholders for the short and long root branches
cropped_tree_str <- gsub(pattern = short, replacement = "short", x = cropped_tree_str)
cropped_tree_str <- gsub(pattern = long, replacement = "long", x = cropped_tree_str)

# iterate and replace each combination of long and short root branches
for (i in brlengths) {
    cropped_tree_str_i <- cropped_tree_str
    cropped_tree_str_i <- gsub(pattern = "short", replacement = i, x = cropped_tree_str_i)
    cropped_tree_str_i <- gsub(pattern = "long", replacement = i + diff_brlength, x = cropped_tree_str_i)
    writeLines(text = cropped_tree_str_i, con = paste("cropped_tree", i, ".newick", sep = ""))
}
