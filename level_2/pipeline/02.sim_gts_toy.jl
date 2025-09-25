"""
simulate_gts_calc_cf_startingtree.jl 
Script for simulating gene trees under the coalescent
using a network in extnewick (input), calculating the CF
table and returning it (output).

Example:
 julia simulate_gts_calc_cf_startingtree.jl
"""

n_reps=30

using Pkg
Pkg.activate("./")

using PhyloNetworks
using PhyloCoalSimulations
using Random
using CSV
using DataFrames
using Base.Threads


### Simulate gene trees
pars = CSV.read("../data/pars.csv",DataFrame)

FILE_LOCK = ReentrantLock()

for par_comb in pars[1,1] # go over all parameter combinations
    par_folder = "../data/pars_$par_comb/"
    par_write_folder= "../data_toy/"
    ngt = 1000

    phys_nos = [1,5,13,15,16,18,25,34]
    Threads.@threads for phy_no in phys_nos

        phy_folder = par_folder*"net_$phy_no/"
        phy_write_folder = par_write_folder*"net_$phy_no/"
        gt_seeds = CSV.read(phy_folder*"seeds.csv",DataFrame)
        input_net = readnewick(phy_folder*"network.extnewick")

        hy_edges = filter(x->x.hybrid==true,input_net.edge)
        PhyloNetworks.setmultiplegammas!(hy_edges,fill(0.5,length(hy_edges)))

        println("we are at setting number $par_comb and phy $phy_no")
        for rep_no in 1:n_reps
            ((rep_no%5) == 0)  && println("rep $rep_no")
            rep_folder = phy_folder*"rep_$rep_no"
            rep_write_folder = phy_write_folder*"rep_$rep_no"
            mkpath(rep_write_folder)
            isfile(rep_write_folder*"/starting_tree.newick") && continue ##we already have a starting tree, move on

            Random.seed!(gt_seeds[rep_no,3])
            # simulate gene trees 
            output_gts = simulatecoalescent(input_net, ngt, 1)

            #write genetrees as .newick format
            writemultinewick(output_gts, rep_write_folder*"/gene_trees.newick")

            ### Calculate concordance factors
            # calculate CFs
            q,t = countquartetsintrees(output_gts;showprogressbar=false)
            df = tablequartetCF(q,t) 
            CSV.write(rep_write_folder*"/CFs.csv", df)

            try
                lock(FILE_LOCK)
                command = `tree-qmc --fast -i $rep_folder/gene_trees.newick -o $rep_folder/starting_tree.newick`
                run(pipeline(command,stdout="out.txt"));
            finally
                unlock(FILE_LOCK)
            end
            rm(rep_write_folder*"/gene_trees.newick") #save memory, delete the gene trees 
            
        end
        GC.gc()
    end
end
rm("out.txt")

