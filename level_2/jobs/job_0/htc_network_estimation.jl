"""
network_estimation_h0.jl

Infer phylogenetic networks using SNaQ

This script takes three arguments, which are in 
their respective order:

-the input directory containing the data files
-the output directory to save results
-the number of threads to use


Usage: 

julia network_estimation.jl input_dir output_dir
"""
# use arguments for feeding tree and cfs
input_loc = ARGS[1]
output_loc = ARGS[2]
nthreads = parse(Int,ARGS[3])

using Distributed

# set up the max number of threads to use based on nruns


# load packages to parallel threads
@everywhere using PhyloNetworks
@everywhere using Functors
@everywhere using DataFrames
@everywhere using CSV


# read the concordance factor table
cfs = readTableCF(input_loc*"CFs.csv")

# read starting tree
global net = readTopology(input_loc*"starting_tree.newick")

#read in parameters for inference
pars = CSV.read(input_loc*"est_pars.csv",DataFrame)

hmax = pars[1,2]
snaq_seed =pars[1,3]
nruns = pars[1,5]

for h in 0:0

snaq!(net, cfs, hmax=h, filename=output_loc*"h_$h", seed=snaq_seed, runs = 1)
global net = readTopology(output_loc*"h_$h.out")

GC.gc()
end

#exit julia
exit()