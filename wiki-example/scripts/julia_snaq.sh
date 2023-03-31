## #!/bin/bash
# input files
file_tree="$1"
file_cf="$2"

#
# extract Julia and julia packages from tar.gz file
tar -xzf julia-1.8.5-linux-x86_64.tar.gz
tar -xzf pro.tar.gz

# add Julia binary to PATH
export PATH=$_CONDOR_SCRATCH_DIR/julia-1.8.5/bin:$PATH

# add Julia packages to DEPOT variable
export JULIA_DEPOT_PATH=$_CONDOR_SCRATCH_DIR/pro

# run Julia script
julia --project=pro julia_snaq.jl "$file_tree" "$file_cf"
