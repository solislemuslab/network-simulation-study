using Pkg
Pkg.activate(".")


using Base.Threads
using PhyloNetworks
using SNaQ
using DataFrames
using CSV




const PRINT_LOCK = ReentrantLock()

function ts_println(args...)
    lock(PRINT_LOCK)
    try
        println(args...)
    finally
        unlock(PRINT_LOCK)
    end
end

############################################
##### Summarize info across replicates #####
############################################
    #Here we are summarizing information that is NOT agnostic of the simulated CFs.
    #So here we will be summarizing info for each replicate of simulated CFs

for par_no in 1:36
    for phy_no in 75:150
        println("on par $par_no phy $phy_no")
        
        net_loc = "../data/pars_$par_no/net_$phy_no/"
        net = readnewick(net_loc*"network.extnewick")
        shrink2cycles!(net)
        shrink3cycles!(net)
        subnet_loc = net_loc*"subnets/"
        mkpath(subnet_loc)

        subnets,subnet_table = PhyloNetworks.getdisplayednetworks(net)
        nrows = size(subnet_table)[1]
        writemultinewick(subnets,subnet_loc*"subnets.newick")
        CSV.write(subnet_loc*"subnet_table.csv",subnet_table)

        subnet_levels = repeat([-1],nrows)
        for (i,subnet) in enumerate(subnets)
            subnet_levels[i] = getlevel(subnet)
        end


        #for rep_no in 1:1
        Threads.@threads for rep_no in 1:30
            ts_println("on par $par_no phy $phy_no rep $rep_no")
            liks = repeat([-Inf],nrows)
            liks_optimized = repeat([-Inf],nrows)
            cf_loc = net_loc*"rep_$rep_no/"
            cfs = readtableCF(cf_loc*"CFs.csv")

            for (i,(subnet,level)) in enumerate(zip(subnets,subnet_levels)) ##look at all subnets
                #println("were at subnet $i")
                if level<2
                    try
                        liks[i] = topologyQpseudolik!(subnet,cfs)
                        liks_optimized[i] = topologymaxQpseudolik!(subnet,cfs).fscore
                    catch
                        liks[i]= -1.0
                        liks_optimized[i] = -1.0
                    end
                else
                    liks[i]= -1.0
                    liks_optimized[i] = -1.0
                end
            end

            subnet_liks = deepcopy(subnet_table)
            insertcols!(subnet_liks, :liks =>liks, :liks_opt => liks_optimized)

            subnet_liks[!,:liks]=liks
            subnet_liks[!,:liks_opt]=liks_optimized
            CSV.write(cf_loc*"liks.csv",subnet_liks)

        end

    end


end