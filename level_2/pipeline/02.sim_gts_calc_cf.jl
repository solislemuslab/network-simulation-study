"""
simulate_gts_calc_cf_startingtree.jl 
Script for simulating gene trees under the coalescent
using a network in extnewick (input), calculating the CF
table and returning it (output).

Example:
 julia simulate_gts_calc_cf_startingtree.jl
"""

n_phy = parse(Int,ARGS[1])
n_reps = parse(Int,ARGS[2])

using Pkg
Pkg.activate("./")

using PhyloNetworks
using PhyloCoalSimulations
using Random
using CSV
using DataFrames


### Simulate gene trees
pars = CSV.read("../data/pars.csv",DataFrame)

for par_comb in pars[1:end,1] # go over all parameter combinations
    par_folder = "../data/pars_$par_comb/"
    ngt = pars[par_comb,5]

    for phy_no in 1:n_phy
        phy_folder = par_folder*"net_$phy_no/"
        gt_seeds = CSV.read(phy_folder*"seeds.csv",DataFrame)
        input_net = readTopology(phy_folder*"network.extnewick")

        println("we are at setting number $par_comb and phy $phy_no")
        for rep_no in 1:n_reps
            ((rep_no%5) == 0)  && println("rep $rep_no")
            rep_folder = phy_folder*"rep_$rep_no/"
            mkpath(rep_folder)

            Random.seed!(gt_seeds[rep_no,3])
            # simulate gene trees 
            output_gts = simulatecoalescent(input_net, ngt, 1)

            #write genetrees as .newick format
            writeMultiTopology(output_gts, rep_folder*"gene_trees.newick")

            ### Calculate concordance factors
            # calculate CFs
            q,t = countquartetsintrees(output_gts;showprogressbar=false)
            df = writeTableCF(q,t) 
            CSV.write(rep_folder*"CFs.csv", df)
            
        end
        GC.gc()
    end
end



