using Pkg
Pkg.activate("./")

using PhyloNetworks
using PhyloCoalSimulations
using QuartetNetworkGoodnessFit
using Random
using CSV
using DataFrames
using Distributed
using Statistics
using SNaQ

nprocs = parse(Int,ARGS[1])

par_file = "seed_pars.csv"
pars = collect(CSV.read(par_file,DataFrame)[1,:])

addprocs(nprocs)

@everywhere using PhyloNetworks, PhyloCoalSimulations, QuartetNetworkGoodnessFit, Random

main_rng = MersenneTwister(pars[4])
ngt = pars[7]
all_seeds = rand(main_rng, UInt64, ngt) 

net = readnewick("network.extnewick")

@everywhere function run_simulation(seed, net)
    Random.seed!(seed)
    return simulatecoalescent(net, 1, 1)[1]
end

output_gts = pmap(seed -> run_simulation(seed, net), all_seeds)

writemultinewick(output_gts, "gene_trees.newick")
command = `tree-qmc --fast -i gene_trees.newick -o qmc_tree.newick`
run(pipeline(command,stdout=devnull));
rm("gene_trees.newick") #save memory, delete the gene trees


q,t = countquartetsintrees(output_gts;showprogressbar=false)
df = DataFrame(tablequartetCF(q,t))
CSV.write("CFs.csv", df)
df=readtableCF(df)

output_gts = nothing #free memory
@everywhere GC.gc()
exit()

#=


##Get a sense of the 'quality of the data' - how well CFs fit the true network
#true_res = quarnetGoFtest!(net,df,false;nsim=10);
#sim_zvals = sort!(true_res[6])
#pval = mean(sim_zvals .>= true_res[2])
#CSV.write("out/true_GoF.csv", DataFrame(pval=pval,zscore=true_res[2],sigma=true_res[3]))

function lesstwo(net)
	return getlevel(net) < 2
end

##Fit CFs to the subnetworks of the true network
subnets,subnet_table = PhyloNetworks.getdisplayednetworks(net;restriction=lesstwo)
nrows = size(subnet_table)[1]


@everywhere using SNaQ
@everywhere function process_subnet((i, subnet), df)   
    level = getlevel(subnet)
        raw_lik = -1.0
        opt_lik = -1.0
        if level < 2
                try
                        raw_lik = topologyQpseudolik!(subnet, df)
                        opt_lik = topologymaxQpseudolik!(subnet, df).fscore
                catch
                end
        end
    return [raw_lik, opt_lik,level,i]   
end
combined_results = pmap(((i, subnet),) -> process_subnet((i, subnet), df), enumerate(subnets))                                                                                                                                                                      
  
                                                                                                                                                                     
column_names = ["raw_lik", "opt_lik", "level", "index"]                                                                                                                    
results_matrix = DataFrame(transpose(hcat(combined_results...)),column_names) 

res=hcat(subnet_table, results_matrix)
@everywhere GC.gc()
CSV.write("out/liks.csv",res)

filter_cond = (res.level .< 2) .& (res.opt_lik .!= -1) 
filtered_df = res[filter_cond, :]
min_opt_lik_row_idx = Int(filtered_df.index[argmin(filtered_df.opt_lik)])

best_subnet = subnets[min_opt_lik_row_idx]
vals=[-1.0,-1.0,-1.0]
try #strange qnet error 
        #TODO look into
        #e.g. par1 phy1 rep 1 h4
        #seems a be an issue with the topology.
        est_res = quarnetGoFtest!(best_subnet,df,false;nsim=1000);
        sim_zvals = sort!(est_res[6])
        pval = mean(sim_zvals .>= est_res[2])
        vals=[pval,est_res[2],est_res[3]]
catch

end

gof_subnet = DataFrame(
    "p_val" => [vals[1]],
    "zscore" => [vals[2]],
    "sigma" => [vals[3]])

CSV.write("out/GoF_best_subnet.csv",gof_subnet)

=#
