using Distributed
addprocs(4)
@everywhere using PhyloNetworks

using PhyloPlots
using CSV

# load input trees
trees = readMultiTopology("hybridLambda_output/net_h1_n6_ho_coal_unit")
# terminals must be renamed in order to match gene-tree names from hybridLambda
# also, an arbitrary position for the hybrid-induced terminal was chosen to be sister to G
start_tree = readTopology("((((A_1,B_1),(C_1,D_1)),(E_1,(F_1,(G_1,Hhyb0.5_1)))),H_1);")

# calculate CFs
q,t = countquartetsintrees(trees)
df = writeTableCF(q,t)
cfs = readTableCF(df)

# calculate a h=0 network
net0 = snaq!(start_tree, cfs, hmax=0, filename="h1_n6_ho_net0", seed=1234)
