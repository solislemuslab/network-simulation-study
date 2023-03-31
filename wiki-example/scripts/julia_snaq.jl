input_net_path = ARGS

using Distributed
addprocs(4, exeflags="--project=$(Base.active_project())")

@everywhere using PhyloNetworks
@everywhere using Functors
@everywhere using DataFrames
@everywhere using CSV


# read the input netowrk in extnewick
gt_data = input_net_path[1]
genetrees = readMultiTopology(gt_data)
tre = genetrees[3]

# Concordance factors from obtained from the gene trees
cf1 = input_net_path[2]
raxmlCF = readTableCF(cf1)

# Snaq
net0 = snaq!(tre,  raxmlCF, hmax=0, filename="net0_bucky", seed=123,runs=4)
net1 = snaq!(net0, raxmlCF, hmax=1, filename="net1_bucky", seed=456,runs=4)

print(net1)# final output that will stored in julread.out
