using Distributed
addprocs(15)
@everywhere using PhyloNetworks
#using PhyloNetworks
using PhyloPlots
using CSV
using DataFrames

gtfilename = ARGS[1]

println("Processing $gtfilename...")

# load input trees
trees = readMultiTopology(gtfilename)

# calculate CFs
q,t = countquartetsintrees(trees)
df = writeTableCF(q,t)
cfs = readTableCF(df)

# add dummy columns to df so that it is get-pop-tree.pl-compliant
zeros = repeat([0], nrow(df))

df_zeros = DataFrame(CF12_34_lo=zeros,
                     CF12_34_hi=zeros,
                     CF13_24_lo=zeros,
                     CF13_24_hi=zeros,
                     CF14_23_lo=zeros,
                     CF14_23_hi=zeros)

df_ticr = hcat(select(df, 1:5),
               select(df_zeros, 1:2),
               select(df, 6),
               select(df_zeros, 3:4),
               select(df, 7),
               select(df_zeros, 5:6))

CSV.write("CFs_$gtfilename.csv", df_ticr)

# calculate the starting tree using part of the TICR pipeline
run(Cmd(["/home/balleng/programs/TICR/scripts/./get-pop-tree.pl", "CFs_$gtfilename.csv"]))

# read starting tree
start_tree = readTopology("CFs_$gtfilename.csv.QMC.tre")

# calculate a h=0 network
net0 = snaq!(start_tree, cfs, hmax=0, filename="snaq_h0_$gtfilename", seed=1234, runs = 15)

# calculate a h=1 network
net1 = snaq!(net0, cfs, hmax=1, filename="snaq_h1_$gtfilename", seed=1234, runs = 15)

# calculate a h=2 network
net2 = snaq!(net1, cfs, hmax=2, filename="snaq_h2_$gtfilename", seed=1234, runs = 15)

#exit julia
exit()
