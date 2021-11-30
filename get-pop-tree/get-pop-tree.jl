#!/usr/bin/julia
using DataFrames
using CSV

"""
This is what get-pop-tree does:

1. Read the qmc input file:
1.1. Get rid of headers
1.2. Use ',' to separate fields
1.3. Create four fields, one for each taxon in the first four columns of the input file, then do the same with the CF frequs in columns 4, 7, and 10 (i.e., exclude the CIs)
1.4. Pick each line which has both quartet configuration and the relative frequencies of ab|cd, ac|bd, and ad|cb.
1.5 Pick the one with the largest value, or both of there are two equally large, or three of all of these are equally well represented f_q = 1/3 for each q, and append it to the qmc input file sequentially and just separeted from the one before by a space.
1.6. Now we have the same information in the format which QMC likes, that is, an array of quartets separated from one another by spaces, in the notation a,b|c,d 
2. Now, QMC seems to need ids instead of actual taxon names, so we need to replace them in the input file before running QMC.
3. Run the program with the modified input and an arbitrary output name which will contain the tree
4. Open the tree file and replace the ids with the actual taxon names
"""

function get_pop_tree(quartet_file, columns, qmc_bin)
    quartets = CSV.read(quartet_file, DataFrame)
end

# simulate arguments below

quartet_file = "1_seqgen.CFs.csv"

columns =  [1, 2, 3, 4, 5, 8, 11]

qmc_bin = "find-cut-Mac"


# read the CFs as a data frame
quartets = CSV.read(quartet_file, DataFrame)

# pick taxa from the first four columns in the data frame and select only unique values
taxa = unique(vcat(quartets[!, columns[1]],
                   quartets[!, columns[2]],
                   quartets[!, columns[3]],
                   quartets[!, columns[4]]))

# create a dict for taxon vs. integer id
dict_taxa = Dict(taxa .=> 1:length(taxa))

# prepare an array for in 
max_cfs = String[]

# collect the largest of the CFs, include more than one if there are ties
for i in eachrow(quartets)
#    println(i[[5, 8, 11]])
    to_choose = maximum(i[[5, 8, 11]]) .== collect(i[[5, 8, 11]])
#    println(to_choose)
    for j in names(i[[5, 8, 11]][to_choose])
#        println(j)
        j = j[[3,4,6,7]]
        println(j)
        if j == "1234"
            println(string(i.taxon1, ",", i.taxon2, "|", i.taxon3, ",", i.taxon4, "\n"))
            push!(max_cfs, string(i.taxon1, ",", i.taxon2, "|", i.taxon3, ",", i.taxon4))
        elseif j == "1324"
            println(string(i.taxon1, ",", i.taxon3, "|", i.taxon2, ",", i.taxon4, "\n"))
            push!(max_cfs, string(i.taxon1, ",", i.taxon3, "|", i.taxon2, ",", i.taxon4))
        elseif j == "1423"
            println(string(i.taxon1, ",", i.taxon4, "|", i.taxon2, ",", i.taxon3, "\n"))
            push!(max_cfs, string(i.taxon1, ",", i.taxon4, "|", i.taxon2, ",", i.taxon3))
        else
            println("Unexpected array of taxa")
        end
    end
end

# now join everything as is required by QMC
qmc_input = join(max_cfs, " ")
