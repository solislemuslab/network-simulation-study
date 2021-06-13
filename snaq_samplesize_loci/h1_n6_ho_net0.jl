### h1_n6_ho_net0.jl
using Distributed
addprocs(4)
@everywhere using PhyloNetworks

using PhyloPlots
using CSV

# load input trees
trees = readMultiTopology("hybridLambda_output/net_h1_n6_ho_coal_unit")
# terminals must be renamed in order to match gene-tree names from hybridLambda
# also, an arbitrary position for the hybrid-induced terminal was chosen to be sister to G
start_tree = readTopology("(15_1,(1_1,((14_1,(((12_1,13_1),(11_1,H3hyb0.5_1)),(7_1,(10_1,(8_1,9_1))))),(((2_1,3_1),(6_1,(5_1,4_1)))))));")

# calculate CFs
q,t = countquartetsintrees(trees)
df = writeTableCF(q,t)
cfs = readTableCF(df)

# calculate a h=0 network
net0 = snaq!(start_tree, cfs, hmax=0, filename="h1_n6_ho_net0", seed=1234)

exit()
