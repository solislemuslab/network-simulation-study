#using Distributed
#addprocs(4)
#@everywhere using PhyloNetworks
using PhyloNetworks
using PhyloPlots
using CSV

gtfilename = ARGS[1]

println("Processing $gtfilename...")

# load input trees
#trees = ARG[2]
trees = readMultiTopology(gtfilename)
#trees = readMultiTopology("gene_trees_cropped_tree_coal_unit")

# calculate CFs
q,t = countquartetsintrees(trees)
df = writeTableCF(q,t)
cfs = readTableCF(df)

CSV.write("CFs_$gtfilename.csv", df)

# calculate the starting tree using part of the TICR pipeline
run(Cmd(["/home/balleng/programs/TICR/scripts/./get-pop-tree.pl", "CFs_$gtfilename.csv"]))

# read starting tree
start_tree = readTopology("CFs_$gtfilename.csv.QMC.tre")

# calculate a h=0 network
net0 = snaq!(start_tree, cfs, hmax=0, filename="snaq_h0_$gtfilename", seed=1234, runs = 1)

# calculate a h=1 network
#net1 = snaq!(net0, cfs, hmax=1, filename="snaq_h1_$gtfilename", seed=1234, runs = 1)

# calculate a h=2 network
#net2 = snaq!(net1, cfs, hmax=2, filename="snaq_h2_$gtfilename", seed=1234, runs = 1)

#exit julia
exit()
