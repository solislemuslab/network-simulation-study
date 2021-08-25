using Distributed
addprocs(4)
@everywhere using PhyloNetworks

using PhyloPlots
using CSV

# load input trees
trees = readMultiTopology("smallnet_ultrametric_coal_unit")
# terminals must be renamed in order to match gene-tree names from hybridLambda
# also, an arbitrary position for the hybrid-induced terminal was chosen to be sister to G
start_tree = readTopology("(E_1:9.583042,(H1hyb_1:5.730323,((A_1:2.353569,B_1:2.353569):4.660416,(C_1:4.035818,D_1:4.035818):2.978167):0.716339):1.852719);")

# calculate CFs
q,t = countquartetsintrees(trees)
df = writeTableCF(q,t)
cfs = readTableCF(df)

# calculate a h=0 network
net0 = snaq!(start_tree, cfs, hmax=0, filename="smallnet_ultrametric_snaq_h0", seed=1234)

# calculate a h=1 network
net1 = snaq!(net0, cfs, hmax=1, filename="smallnet_ultrametric_snaq_h1", seed=1234)

# calculate a h=2 network
net2 = snaq!(net1, cfs, hmax=2, filename="smallnet_ultrametric_snaq_h2", seed=1234)

#exit julia
exit()
