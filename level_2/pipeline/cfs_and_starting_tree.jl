"""
cfs_and_starting_treer.jl

Script for calculating the concordance factor table 
with PhyloNetworks and the initial tree with TICR's
get-pop-tree.pl script which uses QMC.

Usage:

julia cfs_and_starting_tree.jl path-to-genetrees-file.hybridlambda

Note: get-pop-tree.pl must be in the path. It should also be executable
"""

using PhyloNetworks
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
# get-pop-tree.pl must be in the path
run(Cmd(["get-pop-tree.pl", "CFs_$gtfilename.csv"]))

exit()
