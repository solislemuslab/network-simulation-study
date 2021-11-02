using Distributed
addprocs(4)
@everywhere using PhyloNetworks

using PhyloPlots
using CSV

# load input trees
trees = readMultiTopology("gene_trees_cropped_tree_coal_unit")
# terminals must be renamed in order to match gene-tree names from hybridLambda

# calculate CFs
q,t = countquartetsintrees(trees)
df = writeTableCF(q,t)
cfs = readTableCF(df)

CSV.write("CFs.csv", df)

# calculate the starting tree using part of the TICR pipeline
run(Cmd(["/home/balleng/programs/TICR/scripts/./get-pop-tree.pl", "CFs.csv"]))

# read starting tree
start_tree = readTopology("CFs.csv.QMC.tre")

# calculate a h=0 network
net0 = snaq!(start_tree, cfs, hmax=0, filename="smallnet_ultrametric_snaq_h0", seed=1234)

# calculate a h=1 network
net1 = snaq!(net0, cfs, hmax=1, filename="smallnet_ultrametric_snaq_h1", seed=1234)

# calculate a h=2 network
net2 = snaq!(net1, cfs, hmax=2, filename="smallnet_ultrametric_snaq_h2", seed=1234)

#exit julia
exit()
