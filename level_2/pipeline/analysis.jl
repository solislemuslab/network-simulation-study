


using Pkg
Pkg.activate(".")
#Pkg.activate("./data_analysis")
using PhyloNetworks
using SNaQ
using DataFrames
using CSV
using Statistics


####################
### Cluster info ###
####################


for par_no in 1:36
    
    output_dir = "../summarized_results/pars_$par_no/"

    res=DataFrame(phy = Int[],rep=Int[], hmax=Int[],
        hwcd=Int[],
        tp=Int[],tn=Int[],fp=Int[],fn=Int[],
        est_rets=Int[],
        n_broad=Int[],n_narrow=Int[],n_exact=Int[],
        broad_compat=Int[],exact_compat=Int[],
        num_better=Int[], valid_subs=Int[], ##comparing subnetworks to est
        best_hwcd=Int[],best_ll=Float64[])

    found_clus = DataFrame(phy = Int[],rep=Int[], hmax=Int[],
        name=String[],
        broad=Int[],
        narrow=Int[],
        exact=Int[])

    subnet_dat  = DataFrame(phy = Int[],rep=Int[], hmax=Int[],
        subnet_ind  = Int[],
        valid_lik= Bool[],
        delta_ll= Float64[],
        delta_ll_opt=Float64[],
        hwcd=Int[])

    


for phy_no in 1:150
    println("we are at par $par_no and phy $phy_no" ) 
    true_net_loc = "../data/pars_$par_no/net_$phy_no/"
    true_net = readnewick(true_net_loc*"network.extnewick")
    subnets = readmultinewick(true_net_loc*"subnets/subnets.newick")
    



for rep_no in 1:30
    subnet_liks = CSV.read(true_net_loc*"rep_$rep_no/liks.csv",DataFrame)

for h_no in 1:5



est_net_loc = "../output/pars/pars_$par_no/phy_$phy_no/rep_$rep_no/h_$h_no.out"

!isfile(est_net_loc) && continue



est_net = readsnaqnetwork(est_net_loc)

rf_dist = hardwiredclusterdistance(true_net,est_net,false)

taxa = (x -> x.name).(true_net.leaf)

##Compute clusters for hybrid nodes (child edge of hybrid nodes) on true and estimated networks

directedges!(est_net) 
est_clusters = Vector{Vector{Bool}}()
est_hyb_descs = repeat([false],outer=length(taxa))
for hyb_nd in est_net.hybrid
    h_e = getchildedge(hyb_nd) ## get the child edge of the hybrid node 
    clus = hardwiredcluster(h_e,taxa)
    est_hyb_descs = est_hyb_descs .|| clus
    push!(est_clusters,clus)
end

true_clusters = Vector{Vector{Bool}}()
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

found_dat  = DataFrame(phy   = fill(phy_no,nhyb),
                        rep   = fill(rep_no,nhyb),
                        hmax  = fill(h_no,  nhyb),
                        name  = (x->x.name).(true_net.hybrid),
                        broad = Int.(broad_vector),
                        narrow = Int.(narrow_vector),
                        exact = Int.(exact_vector)
                        )
append!(found_clus,found_dat)

##Blob info and level
#=
blobs = PhyloNetworks.blobInfo(true_net)
true_net_level = max((x->length(x)).(blobs[2])...)
true_n_blobs = length(biconnectedComponents(true_net,true))
=#


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


push!(res, (phy_no,rep_no,h_no,
            rf_dist,
            tp,tn,fp,fn,
            length(est_clusters), 
            length(found_clusts), #n_broad
            length(min_matching_clusts),#n_narrow
            sum(length.(exact_clusts).>0),#n_exact
            broad_compat_hybs,
            exact_compat_hybs,
            num_better,
            valid_subs,
            best_rf,
            best_ll
            )
)


end # end h_no
end # end rep_no
end # end phy_no
CSV.write("../summarized_results/pars_$par_no/clusters.csv",res)
CSV.write("../summarized_results/pars_$par_no/found_clades.csv",found_clus)
CSV.write("../summarized_results/pars_$par_no/subnet_dat.csv",subnet_dat)



end



