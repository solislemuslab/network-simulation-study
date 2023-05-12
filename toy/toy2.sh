#!/bin/bash
file_tree="$1"
file_cf="$2"
#
# extract Julia tar.gz file
tar -xzf julia-1.8.5-linux-x86_64.tar.gz
tar -xzf pro.tar.gz

# add Julia binary to PATH
export PATH=$_CONDOR_SCRATCH_DIR/julia-1.8.5/bin:$PATH
#export PATH=$/mnt/ws/home/cacostacortez/julia-1.6.3/bin:$PATH

# add Julia packages to DEPOT variable
export JULIA_DEPOT_PATH=$_CONDOR_SCRATCH_DIR/pro


# run Julia script
# julia --project=pro snaq_julia2_1.jl "$file_tree" "$file_cf"

# run Julia script

# SEED=`cat gt_seed`
SEED=1058

julia --project=pro network_estimation.jl "$file_cf" "$file_tree" 0 4 4 12 > snaq_outgroup_h0.outerr 2>&1
julia --project=pro network_estimation.jl "$file_cf" snaq_output_h0.out 1 4 4 13 > snaq_outgroup_h1.outerr 2>&1
julia --project=pro network_estimation.jl "$file_cf" snaq_output_h1.out 2 4 4 14 > snaq_outgroup_h2.outerr 2>&1
julia --project=pro network_estimation.jl "$file_cf" snaq_output_h2.out 3 4 4 15 > snaq_outgroup_h3.outerr 2>&1
