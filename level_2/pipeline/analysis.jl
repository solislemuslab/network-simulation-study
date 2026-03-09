


using Pkg
Pkg.activate(".")
#using Revise
#Pkg.activate("./data_analysis")
using PhyloNetworks
using SNaQ
using DataFrames
using CSV
using Statistics

include("./aux_functions.jl")


####################
### Cluster info ###
####################

# there are 36 parameter settings (par_no), each with 150 true networks (phy_no), each with 30 replicates (rep_no), each with 5 hmax (h_no) settings
# we loop thru all of these and summarize the results

for par_no in 1:36
    
    
    output_dir = "../summarized_results/pars_$par_no/"

    clus_file = output_dir*"clusters.csv"
    found_clades_file = output_dir*"found_clades.csv"

        # File doesn't exist: Write the HEADER now to avoid race conditions later
        # Create an empty DataFrame with your expected schema
        res=DataFrame(phy = Int[],rep=Int[], hmax=Int[],
            hwcd=Int[],
            tp=Int[],tn=Int[],fp=Int[],fn=Int[], ##hybrid detetction true pos/neg and false pos/neg
            est_rets=Int[], #number of reticulations in the estimated network 
            n_broad=Int[],n_narrow=Int[],n_exact=Int[], #number of (estimated) clusters that map to the true network
            broad_compat=Int[],exact_compat=Int[],
            ave_disp_rf=Float64[],total_disp_rf=Int[], #RF distances for displayed trees 
            gof_pval=Float64[],zscore=Float64[],sigma = Float64[],
            CF_dist=Float64[],
            canon_hwcd=Int[],fu_hwcd=Int[],
            #quar_found = Float64[], ##proportion of true quarnets found in est net
            #quar_compat = Float64[] ##proportion of est quarnets compatible with the true net
            #num_better=Int[], valid_subs=Int[], ##comparing subnetworks to est
            #best_hwcd=Int[],best_ll=Float64[]
        )
        found_clus = DataFrame(phy = Int[],rep=Int[], hmax=Int[], 
        name=String[], #name of the hybrid node in the true network
        broad=Int[], #number of clusters compatible with this hybrid was found. Compatible in that an est cluster has a subset of the true cluster
        narrow=Int[], #the number of clusters that map to this true cluster with minimum distance. Distance being the difference in the number of taxa between the two clusters 
        exact=Int[]) #the number of clusters that map exactly to this true cluster


    if isfile(clus_file)
        # File exists: Load done IDs
        ids_df = CSV.read(clus_file, DataFrame; select=[:phy, :rep, :hmax])
        processed_ids = Set([(r.phy, r.rep, r.hmax) for r in eachrow(ids_df)])
    else
        #create empty files
        CSV.write(clus_file, res,writeheader=true)
        CSV.write(found_clades_file, found_clus,writeheader=true)
        processed_ids = Set{Tuple{Int,Int,Int}}()
    end


  
    #=
    subnet_dat  = DataFrame(phy = Int[],rep=Int[], hmax=Int[],
        subnet_ind  = Int[],
        valid_lik= Bool[],
        delta_ll= Float64[],
        delta_ll_opt=Float64[],
        hwcd=Int[])
    =#
    

df_lock = ReentrantLock()
Threads.@threads for phy_no in 1:150
#for phy_no in 1:150

    phy_res = deepcopy(res)
    phy_found_clus = deepcopy(found_clus)


    println("we are at par $par_no and phy $phy_no" ) 
    true_net_loc = "../data/pars_$par_no/net_$phy_no/"
    true_net = readnewick(true_net_loc*"network.extnewick")
    true_taxa = sort!(tiplabels(true_net))
    true_m = hardwiredclusters(true_net,true_taxa)
    #subnets = readmultinewick(true_net_loc*"subnets/subnets.newick")

    
    #canonical network
    canon_net = PhyloNetworks.canonicalnetwork(true_net)
    canon_m = hardwiredclusters(canon_net,true_taxa)
    
    
    #get FU stable network
    true_mtree = PhyloNetworks.unfold_network(true_net)
    fu_net = PhyloNetworks.stablefold_multree(true_mtree)
    level_set, name_set = PhyloNetworks.getleveltraversal(fu_net)
    true_fold_m = hardwiredclusters(fu_net,true_taxa)
    

    
    #= NO SQUIRREL THINGS WHEN MULTITHREADING 
    #quarnets 
    sq_net = newick_to_squirrel_semi_network(true_net)
    sq_qnets = sq_net.quarnets()
    =#

    #only consider displayed trees if the number of hybrids is small enough
    do_displayed = false
    if true_net.numhybrids <= 10
        do_displayed = true
        true_displayed_trees = PhyloNetworks.getdisplayednetworks(true_net;restriction=only_trees)[1]
        true_displayed_m = (x-> hardwiredclusters(x,true_taxa)).(true_displayed_trees)
    end

    #parent_trees = PhyloNetworks.getweaklydisplayedtrees(true_net)

    


for rep_no in 1:30
    #subnet_liks = CSV.read(true_net_loc*"rep_$rep_no/liks.csv",DataFrame)
    !isfile("../output/pars/pars_$par_no/phy_$phy_no/rep_$rep_no/CFs.csv") && continue

    # Redirect to devnull
    #cfs = redirect_stdout(devnull) do
    #    readtableCF("../output/pars/pars_$par_no/phy_$phy_no/rep_$rep_no/CFs.csv")
    #end
    cfs = readtableCF("../output/pars/pars_$par_no/phy_$phy_no/rep_$rep_no/CFs.csv")

for h_no in 1:5


##skip if already processed
if (phy_no,rep_no,h_no) ∈ processed_ids
    #println("Skipping $rep_no hmax $h_no as already processed")
    continue
end
#println("rep $rep_no hmax $h_no")
##create empty dataframes to store results


est_net_loc = "../output/pars/pars_$par_no/phy_$phy_no/rep_$rep_no/h_$h_no.out"

!isfile(est_net_loc) && continue

#true_net_loc = "../output/pars/pars_$par_no/phy_$phy_no/rep_$rep_no/"
#true_net = readnewick(true_net_loc*"network.extnewick")

gof_dat = CSV.read("../output/pars/pars_$par_no/phy_$phy_no/rep_$rep_no/GoF_$h_no",DataFrame)

est_net = readsnaqnetwork(est_net_loc)
CF_dist = -3.0 #this should never happen
try
    topologyQpseudolik!(est_net, cfs)
    df_wide = fittedquartetCF(cfs)


    CF_dist= 
    mean(abs.(df_wide[:,5] .- df_wide[:,8]).+ ##obs12 - exp12 
        abs.(df_wide[:,6] .- df_wide[:,9]).+ ##obs13 - exp13
        abs.(df_wide[:,7] .- df_wide[:,10])) ##obs14 - exp14
catch e
    println("Error computing expected CFs for par $par_no phy $phy_no rep $rep_no hmax $h_no")
    #println(e)
    CF_dist = -1.0
end


    ######compute RF distance between true and estimated displayed only_trees
    if do_displayed
        if est_net.numhybrids == 0
            est_displayed_trees = [est_net]
        else
            est_displayed_trees = PhyloNetworks.getdisplayednetworks(est_net;restriction = only_trees)[1]
        end
        min_dists = Int[]
        for (true_disp,true_disp_m) in zip(true_displayed_trees,true_displayed_m)
            dist_min = Inf
            for est_disp in est_displayed_trees
                d = hardwiredclusterdistance_firstrooted!(true_disp,est_disp,true_taxa;M1=true_disp_m)
                if d < dist_min
                    dist_min = d
                end
            end
            push!(min_dists,dist_min)
        end
        ave_disp_rf = mean(min_dists)
        total_disp_rf = sum(min_dists)

    else
        ave_disp_rf = -1.0
        total_disp_rf = -1.0
    end



###compute RF distance between reduced versions of the networks
canon_est = PhyloNetworks.canonicalnetwork(est_net)
canon_hwcd = hardwiredclusterdistance_firstrooted!(canon_net,
canon_est,
true_taxa;
M1=canon_m)

if est_net.numhybrids==0
    fu_est = est_net
else
    est_mtree = PhyloNetworks.unfold_network(est_net)
    fu_est = PhyloNetworks.stablefold_multree(est_mtree)
end
fu_hwcd = hardwiredclusterdistance_firstrooted!(fu_net,
fu_est,
true_taxa;
M1=true_fold_m)

#squirrel results can't be multithreaded safely - It's a Python thing
#=
sq_est = newick_to_squirrel_semi_network(est_net)
sq_est_qnets = sq_est.quarnets()
quars_found = consistency_score(sq_qnets,sq_est_qnets)
quars_consistent = consistency_score(sq_est_qnets,sq_qnets)
=#





#=
#get possible foldings of the unfolded network
mtree = readnewick(writenewick(PhyloNetworks.unfold_network(est_net)))
level_set, name_set = PhyloNetworks.getleveltraversal(mtree)

index_map = PhyloNetworks.getinextendible_nodes(level_set,name_set)
updown_links = PhyloNetworks.link_inextnodes(name_set,index_map)
fo = PhyloNetworks.foldorders(index_map,updown_links)
#the FU network is the fo with the fewest foldings
fo_prime = fo[argmin(length.(fo))]

fold_net = PhyloNetworks.fold_multree(mtree,
        fo_prime,
        updown_links,
        index_map)
est_fold_net = PhyloNetworks.canonicalnetwork!(readnewick(writenewick(fold_net)))##Hacky fix for badly confirmed network 
FU_hwcd =hardwiredclusterdistance_firstrooted!(true_fold_net,est_fold_net,true_taxa;M1=true_fold_m)
################
=#

if sort(tipLabels(true_net)) != sort(tipLabels(est_net))
    println("Taxa mismatch at rep $rep_no hmax $h_no, skipping")
    continue
end

rf_dist = hardwiredclusterdistance_firstrooted!(true_net,est_net,true_taxa;M1=true_m)

taxa = (x -> x.name).(true_net.leaf)

##Compute clusters for hybrid nodes (child edge of hybrid nodes) on true and estimated networks

directedges!(est_net) 
est_clusters = Vector{Vector{Bool}}() # the first vector iterates over hybrid nodes, the second denotes taxa that are present in the cluster
est_hyb_descs = repeat([false],outer=length(taxa))
for hyb_nd in est_net.hybrid
    h_e = getchildedge(hyb_nd) ## get the child edge of the hybrid node 
    clus = hardwiredcluster(h_e,taxa)
    est_hyb_descs = est_hyb_descs .|| clus
    push!(est_clusters,clus)
end

true_clusters = Vector{Vector{Bool}}() # the first vector iterates over hybrid nodes, the second denotes taxa that are present in the cluster
true_hyb_descs = repeat([false],outer=length(taxa))
true_hyb_names = String[]
for hyb_nd in true_net.hybrid
    h_e = getchildedge(hyb_nd)
    clus = hardwiredcluster(h_e,taxa)
    true_hyb_descs = true_hyb_descs .|| clus
    push!(true_clusters, clus)
    push!(true_hyb_names,hyb_nd.name)
end
nhyb = length(true_net.hybrid)

##Find rates of (mis)identification of hybrid descent
tp = sum(true_hyb_descs .&& est_hyb_descs)
fp = sum(.!true_hyb_descs .&& est_hyb_descs)
tn = sum(.!true_hyb_descs .&& .!est_hyb_descs)
fn =  sum(true_hyb_descs .&& .!est_hyb_descs)

##find compatible true clusters for each estimated cluster
compat_clusts = Vector{Vector{Int}}()
dist_clusts = Vector{Vector{Int}}() #'distance' from the true cluster to the estimated one 
for est_cluster in est_clusters
    compat_with_est = Vector{Int}() #Index 
    dist_with_est = Vector{Int}()
    for true_ind in 1:length(true_clusters)
        true_cluster = true_clusters[true_ind]
        if all(true_cluster .|| (.!true_cluster .&& .!est_cluster)) ## The estimated cluster doesn't have any tips that the true doesn't
            push!(compat_with_est,true_ind)
            clust_dist = sum(true_cluster .&& .!est_cluster) ## how many tips do they not have in common
            push!(dist_with_est,clust_dist)
        end
    end
    push!(compat_clusts,compat_with_est)
    push!(dist_clusts,dist_with_est)
end

found_clusts = unique(vcat(compat_clusts...)) ## gracious estimate of found clusters
min_matching_clusts = filter(x -> !isempty(x[2]), collect(zip(compat_clusts,dist_clusts)))
min_matching_clusts = (x -> x[1][findmin(x[2])[2]]).(min_matching_clusts) #conservative but not nessecarily accurate mapped found clusters.  
exact_found_clusts = ((x,y) -> x[y.==0]).(compat_clusts,dist_clusts)
exact_clusts = unique(vcat(exact_found_clusts...))

#Number of hybs in est_net that map (in the broad/exact sense) to a true hyb_nd
broad_compat_hybs = sum(.!isempty.(compat_clusts))
exact_compat_hybs = sum(.!isempty.(exact_found_clusts)) 

broad_vector = zeros(nhyb)
broad_vector[found_clusts] .= 1

narrow_vector = zeros(nhyb)
narrow_vector[min_matching_clusts] .= 1

exact_vector = zeros(nhyb)
exact_vector[exact_clusts] .=1


##Blob info and level
#=
blobs = PhyloNetworks.blobInfo(true_net)
true_net_level = max((x->length(x)).(blobs[2])...)
true_n_blobs = length(biconnectedComponents(true_net,true))
=#

#=
#loglik and RF differences from subnetworks
n_subnets = length(subnets)
rfs =repeat([-1],n_subnets)
delta_ll = repeat([-Inf],n_subnets)
delta_ll_opt = repeat([-Inf],n_subnets)
valids = repeat([false],n_subnets) 
for (i,subnet) in enumerate(subnets)
    rfs[i] = hardwiredclusterdistance(subnet,est_net,false)
    valids[i] = subnet_liks.liks[i] >= 0  #check if a valid liklihood
    if valids[i] 
        delta_ll[i] = est_net.fscore - subnet_liks.liks[i]
        delta_ll_opt[i] = est_net.fscore - subnet_liks.liks_opt[i]
    else
        delta_ll[i] = -Inf
        delta_ll_opt[i] = -Inf
    end
end

subnet_deltas  = DataFrame(phy   = fill(phy_no,n_subnets),
                        rep   = fill(rep_no,n_subnets),
                        hmax  = fill(h_no,  n_subnets),
                        subnet_ind  = 1:n_subnets,
                        valid_lik= valids,
                        delta_ll=delta_ll,
                        delta_ll_opt=delta_ll_opt,
                        hwcd=rfs
                        )
append!(subnet_dat,subnet_deltas)


best_ind = argmax(delta_ll_opt)
best_ll = delta_ll_opt[best_ind]
best_rf = rfs[best_ind]
num_better = sum(delta_ll_opt[valids] .< 0.0 )
valid_subs = sum(valids)

=#

    push!(phy_res, (phy_no,rep_no,h_no,
                rf_dist,
                tp,tn,fp,fn,
                length(est_clusters), 
                length(found_clusts), #n_broad
                length(min_matching_clusts),#n_narrow
                length(exact_clusts),#n_exact
                broad_compat_hybs,
                exact_compat_hybs,
                ave_disp_rf,
                total_disp_rf,
                gof_dat[1,1], #pval
                gof_dat[1,2], #zscore
                gof_dat[1,3], #sigma
                CF_dist,
                canon_hwcd,
                fu_hwcd,
                #quars_found,
                #quars_consistent,
                # num_better,
                # valid_subs,
                # best_rf,
                # best_ll
                )
    )

    found_dat  = DataFrame(phy   = fill(phy_no,nhyb),
                        rep   = fill(rep_no,nhyb),
                        hmax  = fill(h_no,  nhyb),
                        name  = (x->x.name).(true_net.hybrid),
                        broad = Int.(broad_vector),
                        narrow = Int.(narrow_vector),
                        exact = Int.(exact_vector)
                        )
    append!(phy_found_clus,found_dat)


end # end h_no
end # end rep_no

    lock(df_lock) do
        CSV.write(clus_file, phy_res; append=true, writeheader=false)
        CSV.write(found_clades_file, phy_found_clus; append=true, writeheader=false)
    end 
end # end phy_no
#CSV.write("../summarized_results/pars_$par_no/clusters.csv",res)
#CSV.write("../summarized_results/pars_$par_no/found_clades.csv",found_clus)
#CSV.write("../summarized_results/pars_$par_no/subnet_dat.csv",subnet_dat)



end

exit()

#########################
####Squirrel analyses####
#########################


for par_no in 1:36
    
    output_dir = "../summarized_results/pars_$par_no/"

    sq_file = output_dir*"squirrel.csv"
    
    sq_res=DataFrame(phy = Int[],rep=Int[], hmax=Int[],
        quar_found = Float64[], ##proportion of true quarnets found in est net
        quar_compat = Float64[] ##proportion of est quarnets compatible with the true net
    )
    
    if isfile(sq_file)
        # File exists: Load done IDs
        ids_df = CSV.read(sq_file, DataFrame; select=[:phy, :rep, :hmax])
        processed_ids = Set([(r.phy, r.rep, r.hmax) for r in eachrow(ids_df)])
    else
        #create empty files
        CSV.write(sq_file, sq_res,writeheader=true)
        processed_ids = Set{Tuple{Int,Int,Int}}()
    end

for phy_no in 1:150

    phy_sq_res = deepcopy(sq_res)

    println("we are at par $par_no and phy $phy_no" ) 
    true_net_loc = "../data/pars_$par_no/net_$phy_no/"
    true_net = readnewick(true_net_loc*"network.extnewick")


    #quarnets
    valid_quars = true 
    sq_net = nothing
    sq_qnets = nothing
    try
        sq_net = newick_to_squirrel_semi_network(true_net)
        sq_qnets = sq_net.quarnets()
    catch e
        println("Error computing quarnets for par $par_no phy $phy_no")
        valid_quars = false
    end

for rep_no in 1:30
    !isfile("../output/pars/pars_$par_no/phy_$phy_no/rep_$rep_no/CFs.csv") && continue

for h_no in 1:5


##skip if already processed
if (phy_no,rep_no,h_no) ∈ processed_ids
    #println("Skipping $rep_no hmax $h_no as already processed")
    continue
end

est_net_loc = "../output/pars/pars_$par_no/phy_$phy_no/rep_$rep_no/h_$h_no.out"

!isfile(est_net_loc) && continue


est_net = readsnaqnetwork(est_net_loc)

##populate results if we can't compute quarnets
if !valid_quars
    push!(phy_sq_res, (phy_no,rep_no,h_no,
            -1.0,
            -1.0
            )
    )
    continue
end

sq_est = newick_to_squirrel_semi_network(est_net)
sq_est_qnets = sq_est.quarnets()
quars_found = consistency_score(sq_qnets,sq_est_qnets)
quars_consistent = consistency_score(sq_est_qnets,sq_qnets)

push!(phy_sq_res, (phy_no,rep_no,h_no,
            quars_found,
            quars_consistent
            )
)

end # end h_no
end # end rep_no
    CSV.write(sq_file, phy_sq_res; append=true, writeheader=false)
end # end phy_no


end
