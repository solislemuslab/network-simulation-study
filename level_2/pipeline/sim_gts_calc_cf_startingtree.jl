"""
simulate_gts.jl 
Script for simulating gene trees under the coalescent
using a network in extnewick (input) and returning the
gene trees in newick (output).
This script reads the first argument as the path to 
the input file, and the second as the path to the 
output file. Arguments nsims and seed are integers
indicating the number of gene trees to simulate and
the random seed respectively.
example
 julia simulate_gts.jl input.extnewick output.newick nsims seed
"""

directory = ARGS[1]
nsims = parse(Int64, ARGS[2])[1]

using PhyloNetworks
using PhyloCoalSimulations
using Random
using CSV
using DataFrames

cd(directory)

network_dirs = readdir()

for i in network_dirs
    # visit each dir containing one network
    cd(i)

    ### Simulate gene trees
    # fetch file names and create necessary variables
    dircontents = readdir()
    network_id = findall( x -> occursin("extnewick", x), dircontents)
    network_file = dircontents[network_id]
    seed_id = findall( x -> occursin("seed", x), dircontents)
    seed_file = dircontents[seed_id]
    output_file = replace(network_file[1], "ext"=>"")
    seed = open(f->read(f, String), "gt_seed")

    # set random seed from the gt_seed file
    Random.seed!(parse(Int64, seed))    

    # read the input netowrk in extnewick
    input_net = readTopology(network_file[1])
    # simulate gene trees 
    output_gts = simulatecoalescent(input_net, nsims, 1)
    # write to file
    writeMultiTopology(output_gts, output_file)

    ### Calculate concordance factors
    # calculate CFs
    q,t = countquartetsintrees(output_gts)
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
    CSV.write("CFs_$output_file.csv", df_ticr)

    ### Use QMC for estimating an initial tree using CFs
    # calculate the starting tree using part of the TICR pipeline
    # get-pop-tree.pl must be in the path
    run(Cmd(["get-pop-tree.pl", "CFs_$output_file.csv"]))

    # go back to data to visit the next directory
    cd("..")
end
