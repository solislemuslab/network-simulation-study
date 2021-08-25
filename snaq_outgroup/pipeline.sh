#!/usr/bin/bash

# simulate a _tree_ in R
Rscript simulate_tree.R

# simulate gene trees on cropped_tree using hybrid-Lambda
hybrid-Lambdav0.6.2 -spcu cropped_tree.newick -num 3000 -seed 1453 -o gene_trees_cropped_tree

# calculate CFs with bucky on the simulated gene trees

# calculateed a tree with QMC on the output from bucky

# run snaq on the gene trees from hybrid-Lambda with the init-tree from QMC 
