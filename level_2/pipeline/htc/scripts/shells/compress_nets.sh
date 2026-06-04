#!/bin/bash

# DAGMan passes the dataset ID
dataset=$1

# CD into the specific job folder
JOB_DIR="/home/justison/jobs/job_${dataset}"
cd $JOB_DIR || exit 1

# SAFETY CHECK 1: If the tarball already exists, we are already done! Skip to exit.
if [ -f "all_nets.tar.gz" ]; then
    echo "Tarball already exists. Skipping compression."
    exit 0
fi

# Create the tarball containing all analysis
tar -czf all_nets.tar.gz h_*.out CFs.csv 

# SAFETY CHECK 2: Only delete the raw files if the tarball was successfully created!
if [ -f "all_nets.tar.gz" ]; then
    rm -f h_*.out
    rm -f CFs.csv
    rm -f qmc_tree.newick
else
    exit 1
fi