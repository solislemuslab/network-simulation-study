


using Pkg
Pkg.activate(".")
#Pkg.activate("./data_analysis")
using PhyloNetworks
using SNaQ
using DataFrames
using CSV
using QuartetNetworkGoodnessFit
using Statistics
using Base.Threads


############################
### Goodness of fit test ###
############################




const PRINT_LOCK = ReentrantLock()
const GOF_LOCK = ReentrantLock()

function ts_println(args...)
    lock(PRINT_LOCK)
    try
        println(args...)
    finally
        unlock(PRINT_LOCK)
    end
end



for par_no in 1:36
    output_dir = "../summarized_results/pars_$par_no/"

    ##Read or make the CSV
    if !isfile(output_dir*"GoF.csv")
        gof_res = DataFrame(phy = Int[],hmax = Int[],pval=Float64[],zscore=Float64[],sigma=Float64[])
    else
        gof_res = CSV.read(output_dir*"GoF.csv",DataFrame)
    end
for phy_no in 1:150
Threads.@threads for rep_no in 1:30
    
        
    for h_no in 1:5
        ts_println("on par $par_no phy $phy_no rep $rep_no h $h_no")

        true_net_loc = "../data/pars_$par_no/net_$phy_no/"
        est_net_loc = "../output/pars/pars_$par_no/phy_$phy_no/rep_$rep_no/h_$h_no.out"
        #check that we are at a valid file location with finished results and that we don't have this file computed
        !isfile(est_net_loc) && continue

        lock(GOF_LOCK)
        try
            if !isempty(gof_res) && any(gof_res.phy .== phy_no) && any(gof_res.hmax .== h_no)
                continue # Skip this iteration if result already exists
            end
        finally
            unlock(GOF_LOCK)
        end


        #read in the estimated and true network
        true_net = readnewick(true_net_loc*"network.extnewick")
        est_net = readnewick(est_net_loc)

        obs_cfs = CSV.read(true_net_loc*"rep_$rep_no/CFs.csv",DataFrame) #CFs from GTs, used to infer the network
        try #strange qnet error 
            #TODO look into
            #e.g. par1 phy1 rep 1 h4
            #seems a be an issue with the topology.
            est_res = quarnetGoFtest!(est_net,obs_cfs,false;nsim=1000);
            sim_zvals = sort!(est_res[6])
            pval = mean(sim_zvals .>= est_res[2])
            vals=[phy_no,h_no,pval,est_res[2],est_res[3]]
        catch
            vals=[phy_no,h_no,-1.0,-1.0,-1.0]
        end


        lock(GOF_LOCK)
        try
            push!(gof_res,)
        finally
            unlock(GOF_LOCK)
        end

        

    end
end #threads loop
    CSV.write(output_dir*"GoF.csv",gof_res)
end
end
