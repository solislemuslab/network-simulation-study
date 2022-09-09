"""
simulate_gts.jl 
Script for simulating gene trees under the coalescent
using a network in extnewick (input) and returning the
gene trees in newick (output).
This script reads the first argument as the path to 
the input file, and the second as the path to the 
output file.
example
 julia extnewick2hybridlambda.jl input.extnewick output.newick
"""

input_net_path = ARGS[1]
output_gts_path = ARGS[2]
nsims = parse(Int64, ARGS[3])
seed = parse(Int64, ARGS[4])

using PhyloNetworks
using PhyloCoalSimulations
using Random

Random.seed!(seed)

# read the input netowrk in extnewick
input_net = readTopology(input_net_path)

# simulate gene trees 
output_gts = simulatecoalescent(input_net, nsims, 1)

# write to file
writeMultiTopology(output_gts, output_gts_path)

exit()
