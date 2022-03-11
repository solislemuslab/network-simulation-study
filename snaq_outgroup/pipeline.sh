#!/usr/bin/bash

# simulate a _tree_ in R
Rscript simulate_tree.R

# simulate gene trees on cropped_tree using hybrid-Lambda
for i in `ls cropped_tree*.newick`; do
    hybrid-Lambdav0.6.2 -spcu $i -num 3000 -seed 1453 -o ${i}_gene_trees
done

# calculate CFs with PhyloNetworks and init_tree with TICR
for i in `ls *.newick_gene_trees_coal_unit`; do
    julia cfs_and_starting_tree.jl $i
done

# prepare the pbs scripts using the template file template.pbs
for i in `ls CFs_*csv`; do
    # create the pbs script from the template file and just append the extension .pbs
    touch $i.pbs
    cat template.pbs > $i.pbs
    # replace the placeholder JOBID with the CF filename
    sed -i "s/JOBID/$i/g" $i.pbs
done

