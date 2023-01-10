## #!/bin/bash
# input files
file_tree="$1"#Gene tree and consensus tree
file_cf="$2"#Concordanse factors
#
# extract Julia and julia packages from tar.gz file
tar -xzf julia-1.6.3-linux-x86_64.tar.gz
tar -xzf my-julia-project.tar.gz

# add Julia binary to PATH
export PATH=$_CONDOR_SCRATCH_DIR/julia-1.6.3/bin:$PATH

# add Julia packages to DEPOT variable
export JULIA_DEPOT_PATH=$_CONDOR_SCRATCH_DIR/my-julia-project

# run Julia script
julia --project=my-julia-project julia_snaq.jl "$file_tree" "$file_cf"
