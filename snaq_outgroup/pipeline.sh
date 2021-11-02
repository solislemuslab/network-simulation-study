#!/usr/bin/bash

# simulate a _tree_ in R
Rscript simulate_tree.R

# simulate gene trees on cropped_tree using hybrid-Lambda
for i in `ls cropped_tree*.newick`; do
    hybrid-Lambdav0.6.2 -spcu $i -num 3000 -seed 1453 -o ${i}_gene_trees
done

# run snaq on the gene trees from hybrid-Lambda with the init-tree from QMC 
julia snaq_cropped_tree.jl
