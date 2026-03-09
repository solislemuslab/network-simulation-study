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

const FILE_LOCK = ReentrantLock()

### Simulate gene trees
pars = CSV.read("../data/pars.csv",DataFrame)

for par_comb in pars[1:end,1] # go over all parameter combinations
    par_folder = "../data/pars_$par_comb/"
    ngt = pars[par_comb,5]

    Threads.@threads for phy_no in 1:n_phy

        phy_folder = par_folder*"net_$phy_no/"
        gt_seeds = CSV.read(phy_folder*"seeds.csv",DataFrame)
        input_net = readnewick(phy_folder*"network.extnewick")

        println("we are at setting number $par_comb and phy $phy_no")
        for rep_no in 1:n_reps
            ((rep_no%5) == 0)  && println("rep $rep_no")
            rep_folder = phy_folder*"rep_$rep_no"
            mkpath(rep_folder)
            #isfile(rep_folder*"/starting_tree.newick") && continue ##we already have a starting tree, move on

            Random.seed!(gt_seeds[rep_no,3])
            # simulate gene trees 
            output_gts = simulatecoalescent(input_net, ngt, 1)

            #write genetrees as .newick format
            writemultinewick(output_gts, rep_folder*"/gene_trees.newick")

            ### Calculate concordance factors
            # calculate CFs
            q,t = countquartetsintrees(output_gts;showprogressbar=false)
            df = tablequartetCF(q,t) 
            CSV.write(rep_folder*"/CFs.csv", df)

            
            try
                lock(FILE_LOCK)
                command = `tree-qmc --fast -i $rep_folder/gene_trees.newick -o $rep_folder/starting_tree.newick`
                run(pipeline(command,stdout="out.txt"));
            finally
                unlock(FILE_LOCK)
            end
            rm(rep_folder*"/gene_trees.newick") #save memory, delete the gene trees 
            
        end
        GC.gc()


    end
end
rm("out.txt")


