
# use arguments for feeding tree and cfs
hmax = parse(Int,ARGS[2])
run_no = parse(Int,ARGS[6])
start_tree = ARGS[4]



using Statistics

# load packages to parallel threads
using PhyloNetworks
using SNaQ
using DataFrames
using CSV
using Random


# read the concordance factor table
cfs = readtableCF("CFs.csv")

# read starting tree
start = readnewick(start_tree)
#read in parameters for inference
pars = collect(CSV.read("seed_pars.csv",DataFrame)[1,:])


snaq_seed =pars[3]
Random.seed!(snaq_seed)
h_seeds = rand(Int64,hmax+1)
Random.seed!(h_seeds[end])
run_seeds= rand(Int64,45+1)
Random.seed!(run_seeds[run_no+1])

snaq!(start,cfs;
        hmax=hmax, filename="h_$(hmax)_run_$(run_no)", seed=run_seeds[run_no+1],
        runs = 1,Nfail=600,probST=0.0) 

#=
GC.gc()
    try #strange qnet error 
        #TODO look into
        #e.g. par1 phy1 rep 1 h4
        #seems a be an issue with the topology.
        est_res = quarnetGoFtest!(net,cfs,false;nsim=1000);
        sim_zvals = sort!(est_res[6])
        pval = mean(sim_zvals .>= est_res[2])
        global vals=[pval,est_res[2],est_res[3]]
    catch
        global vals=[-1.0,-1.0,-1.0]
    end

    gof_subnet = DataFrame(
        "p_val" => [vals[1]],
        "zscore" => [vals[2]],
        "sigma" => [vals[3]])

CSV.write("out/GoF_$h",gof_subnet)


GC.gc()
end
=#
#exit julia
exit()
