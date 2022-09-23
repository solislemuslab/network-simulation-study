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

cd(directory)

network_dirs = readdir(join=true)

for i in network_dirs
    cd(i)
    dircontents = readdir(join=true)
    network_id = findall( x -> occursin("extnewick", x), dircontents)
    network_file = dircontents[network_id]
    seed_id = findall( x -> occursin("seed", x), dircontents)
    seed_file = dircontents[seed_id]
    output_file = replace(network_file[1], "ext"=>"")
    seed = open(f->read(f, String), "gt_seed")

    Random.seed!(parse(Int64, seed))    

    # read the input netowrk in extnewick
    input_net = readTopology(network_file[1])
    # simulate gene trees 
    output_gts = simulatecoalescent(input_net, nsims, 1)
    # write to file
    writeMultiTopology(output_gts, output_file)
    cd("..")
end

exit()
