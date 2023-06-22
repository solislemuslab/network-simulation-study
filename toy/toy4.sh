#!/bin/bash

# extract Julia tar.gz file
tar -xzf julia-1.8.5-linux-x86_64.tar.gz
tar -xzf pro.tar.gz

# add Julia binary to PATH
export PATH=$_CONDOR_SCRATCH_DIR/julia-1.8.5/bin:$PATH

# add Julia packages to DEPOT variable
export JULIA_DEPOT_PATH=$_CONDOR_SCRATCH_DIR/pro

# run Julia script

tar -xzvf files.tar.gz

SEED=`cat gt_seed`

julia --project=pro network_estimation.jl *.csv *.tre 0 4 4 `expr $SEED + 1` > snaq_outgroup_h0.outerr 2>&1
julia --project=pro network_estimation.jl *.csv snaq_output_h0.out 1 4 4 `expr $SEED + 3` > snaq_outgroup_h1.outerr 2>&1
julia --project=pro network_estimation.jl *.csv snaq_output_h1.out 2 4 4 `expr $SEED + 5` > snaq_outgroup_h2.outerr 2>&1
julia --project=pro network_estimation.jl *.csv snaq_output_h2.out 3 4 4 `expr $SEED + 7` > snaq_outgroup_h3.outerr 2>&1
