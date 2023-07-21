"""
network_estimation_h0.jl

Infer phylogenetic networks using SNaQ

This script takes five arguments, which are in
their respective order:

- cf_table: The path to the CSV file containing the CFs
- init_tree: the path to the initial tree, for instance from
  an analysis using maximum likelihood, get-pop-tree.pl,
  or a previous run of SNaQ, in whose .out file we have a tree.
- h: The maximum h in the iteration
- nthreads: The number of threads to use for the analysis. It will
  reserve nthreads+1 for the main Julia process.
- nruns: Number of parallel runs to be carried out
- seed: the random seed

Usage:

julia network_estimation.jl cf_table init_tree h nruns seed
"""

# use arguments for feeding tree and cfs
cf_table = ARGS[1]
init_tree = ARGS[2]
h = parse(Int, ARGS[3])
nthreads = parse(Int, ARGS[4])
nruns = parse(Int, ARGS[5])
seed = parse(Int, ARGS[6])

using Distributed

# set up the max number of threads to use based on nruns
addprocs(nthreads, exeflags="--project=$(Base.active_project())")

# load packages to parallel threads
@everywhere using PhyloNetworks
@everywhere using Functors
@everywhere using DataFrames
@everywhere using CSV

# read the concordance factor table
cfs = readTableCF(cf_table)

# check the size of cfs
println(Base.format_bytes(Base.summarysize(cfs)))

# read starting tree
start_tree = readTopology(init_tree)

# check the size of start_tree
println(Base.format_bytes(Base.summarysize(start_tree)))

# calculate a h=1 network
net = snaq!(start_tree, cfs, hmax=h, filename="snaq_output_h$h", seed=seed, runs = nruns)

# check the size of  net
println(Base.format_bytes(Base.summarysize(net)))

#exit julia
#exit()
