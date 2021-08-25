library(TreeSim)

set.seed(1453)

# simulate a tree
tree <- sim.bd.taxa(n = 30, numbsim = 1, lambda = 1, mu = 0, frac = 1, complete = TRUE, stochsampling = FALSE)

# drop the tips that will create a long branch for the root
drop_group <- c("t30", "t7", "t13", "t15", "t28", "t14", "t11", "t22", "t16", "t26", "t17", "t24", "t4", "t3", "t19", "t21", "t29", "t25")

tree <- plot(drop.tip(tree[[1]], tip = drop_group))
