using PhyloNetworks
using DataFrames
using CSV
using SNaQ
using Statistics
using DataStructures
using QuartetNetworkGoodnessFit
using Distributed

n_workers = parse(Int, ARGS[1])
addprocs(n_workers)

@everywhere using PhyloNetworks
@everywhere using QuartetNetworkGoodnessFit

include("aux_functions.jl")


println("loaded packages")
flush(stdout)

cfs = CSV.read("CFs.csv", DataFrame)
true_net_loc = "network.extnewick"

pars = collect(CSV.read("seed_pars.csv",DataFrame)[1,:])
phy_no= pars[13]
rep_no= pars[14]


###############################
####Preprocess true network####
###############################

true_net = readnewick(true_net_loc)

println("read in core data")

true_taxa = sort!(tiplabels(true_net))
true_m = hardwiredclusters(true_net,true_taxa)

#quarnets
valid_quars = true 
sq_net = nothing
sq_qnets = nothing
try
    global sq_net = newick_to_squirrel_semi_network(true_net)
    global sq_qnets = sq_net.quarnets()
catch e
    println("Error computing quarnets for par $par_no phy $phy_no")
    global valid_quars = false
end


#canonical network
canon_net = PhyloNetworks.canonicalnetwork(true_net)
canon_m = hardwiredclusters(canon_net,true_taxa)


#get FU stable network
true_mtree = PhyloNetworks.unfold_network(true_net)
fu_net = PhyloNetworks.stablefold_multree(true_mtree)
level_set, name_set = PhyloNetworks.getleveltraversal(fu_net)

muledge_merge!(fu_net)

true_fold_m = hardwiredclusters(fu_net,true_taxa)

do_displayed = false
if true_net.numhybrids <= 10
    do_displayed = true
    true_displayed_trees = PhyloNetworks.getdisplayednetworks(true_net;restriction=only_trees)[1]
    true_displayed_m = (x-> hardwiredclusters(x,true_taxa)).(true_displayed_trees)
end




#####################
###make data files###
#####################

clus_file = "clusters.csv"
found_clades_file = "found_clades.csv"
sq_file = "squirrel.csv"

phy_res=DataFrame(phy = Int[],rep=Int[], hmax=Int[],
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
phy_found_clus = DataFrame(phy = Int[],rep=Int[], hmax=Int[], 
    name=String[], #name of the hybrid node in the true network
    broad=Int[], #number of clusters compatible with this hybrid was found. Compatible in that an est cluster has a subset of the true cluster
    narrow=Int[], #the number of clusters that map to this true cluster with minimum distance. Distance being the difference in the number of taxa between the two clusters 
    exact=Int[]) #the number of clusters that map exactly to this true cluster

phy_sq_res=DataFrame(phy = Int[],rep=Int[], hmax=Int[],
    quar_found = Float64[], ##proportion of true quarnets found in est net
    quar_compat = Float64[] ##proportion of est quarnets compatible with the true net
)

println("preprocessed true net")

####################
### Cluster info ###
####################
for h_no in 0:5

    println("on hval: $h_no")
    flush(stdout)

    est_net_loc = "h_$h_no.out"

    est_net = readnewick(est_net_loc)


    try #strange qnet error 
        #TODO look into
        #e.g. par1 phy1 rep 1 h4
        #seems a be an issue with the topology.
        dfcfs = readtableCF("CFs.csv")
        est_res = quarnetGoFtest!(est_net,cfs,false;nsim=1000);
        sim_zvals = sort!(est_res[6])
        pval = mean(sim_zvals .>= est_res[2])
        global vals=[pval,est_res[2],est_res[3]]
    catch
        global vals=[-1.0,-1.0,-1.0]
    end

    gof_dat = DataFrame(
        "p_val" => [vals[1]],
        "zscore" => [vals[2]],
        "sigma" => [vals[3]])


CF_dist = -3.0 #this should never happen
try
    exp_cfs= network_expectedCF(est_net;showprogressbar=false)
    quar_dists=Vector{Float64}(undef,nrow(cfs))
    for rw in 1:nrow(cfs)
        quar_dists[rw] = sum(abs.(collect(exp_cfs[1][rw].data) .- collect(cfs[rw,[6,7,8]])))
    end
    CF_dist=mean(quar_dists)

    #CF_dist= 
    #mean(abs.(df_wide[:,5] .- df_wide[:,8]).+ ##obs12 - exp12 
    #    abs.(df_wide[:,6] .- df_wide[:,9]).+ ##obs13 - exp13
    #    abs.(df_wide[:,7] .- df_wide[:,10])) ##obs14 - exp14
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
    muledge_merge!(fu_est)
end
fu_hwcd = hardwiredclusterdistance_firstrooted!(fu_net,
fu_est,
true_taxa;
M1=true_fold_m)



if sort(tipLabels(true_net)) != sort(tipLabels(est_net))
    println("Taxa mismatch at rep $rep_no hmax $h_no, skipping")
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

####Squirrel things
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



####save results to DF

push!(phy_sq_res, (phy_no,rep_no,h_no,
            quars_found,
            quars_consistent
            )
)


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
end

CSV.write(clus_file, phy_res)
CSV.write(found_clades_file, phy_found_clus)
CSV.write(sq_file, phy_sq_res)

exit()

    

