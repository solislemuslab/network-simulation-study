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

julia network_estimation.jl input_dir output_dir
"""
# use arguments for feeding tree and cfs
input_loc = ARGS[1]
output_loc = ARGS[2]
nthreads = ARGS[3]

using Pkg
Pkg.activate("./")

using Distributed

nthreads = 4



# set up the max number of threads to use based on nruns
addprocs(nthreads)

# load packages to parallel threads
@everywhere using PhyloNetworks
@everywhere using Functors
@everywhere using DataFrames
@everywhere using CSV



for phy in [150,149,148,147,146,145,144,143,142,141,140]
  println(phy)
  input_loc = "../data/pars_31/net_$phy/rep_1/"
  output_loc = "../output/pars_31/net_$phy/rep_1/"

  # read the concordance factor table
  cfs = readTableCF(input_loc*"CFs.csv")

  # read starting tree
  net = readTopology(input_loc*"starting_tree.newick")

  #read in parameters for inference
  pars = CSV.read(input_loc*"est_pars.csv",DataFrame)

  hmax = pars[1,2]
  snaq_seed =pars[1,3]
  nruns = pars[1,5]

  for h in 1:hmax
    net = snaq!(net, cfs, hmax=h, filename=output_loc*"h_$h", seed=snaq_seed, runs = 5)
    GC.gc()
  end

end
#exit julia
#exit()
