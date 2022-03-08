"""
extnewick2hybridlambda.jl 

Script for format conversion between
extended newick (input) and hybrid-Lambda (output)

This script reads the first argument as the path to 
the input file, and the second as the path to the 
output file.

example

$ julia extnewick2hybridlambda.jl input.extnewick input.hybridlambda
"""

input_net_path = ARGS[1]
output_net_path = ARGS[2]

using PhyloNetworks, PhyloPlots
using RCall

# read the input netowrk in extnewick
input_net = readTopology(input_net_path)

# convert to hybridlambda format
output_net = hybridlambdaformat(input_net)

# write to file
open(output_net_path, "w") do io
    write(io, output_net)
end

exit()
