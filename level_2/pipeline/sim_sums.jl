

using Pkg
Pkg.activate(".")
using Base.Threads


using PhyloNetworks
using SNaQ
using DataFrames
using CSV
using Statistics




######################################
##### Summarize info phylogenies #####
######################################
    #Here we are summarizing information about the simulated networks. This is agnostic of the simulated CFs
    #So we will be summarizing info at the phylogeny-level since this is the same regardless of the replicate simualted CF.

for par_no in 1:36
    
    output_dir = "../data/pars_$par_no/"

    hyb_dat= DataFrame(phy=Int[],
    name=String[],blob_no=Int[],
    dist50=Float64[],
    blob_size=Int[],
    blob_level=Int[],
    external = Bool[],
    stacked = Bool[],
    cycle_size=Int[],
    cycle_diam=Float64[],
    cycled=Bool[],
    cy_adjacent=Bool[])
    
    branch_lengths = Float64[]

for phy_no in 1:150
    


    println("we are at par $par_no and phy $phy_no" ) 
    true_net_loc = "../data/pars_$par_no/net_$phy_no/"

    #=
    found_true_net = false
    true_net_loc = ""
    for try_net in 1:30
        trial_net_loc = "../output/pars/pars_$par_no/phy_$phy_no/rep_$try_net/"
        if !isfile(trial_net_loc*"network.extnewick") #move on if the file does not exist
            continue
        else
        true_net_loc = trial_net_loc
        end
    end
    =#
    net = readnewick(true_net_loc*"network.extnewick")
    preorder!(net)
    output_net_dir = output_dir*"net_$phy_no/hybs/"
    mkpath(output_net_dir)

    ##get internal branch lengths
    internal_edges = filter(x-> !getchild(x).leaf,net.edge)
    bl = (x->x.length).(internal_edges)
    bl = filter(x-> !isapprox(x,0),bl) #remove 0 length branches
    append!(branch_lengths,bl)

    hyb_nds = net.hybrid
    hyb_blob_map = Dict{Int,Tuple{Int,Int,Int,Bool,Bool,Int,Float64}}() #connect a blob (cycle size and level) to a node
    hyb_cycle_map= Dict(nd.number => Vector{Tuple{String,Bool}}() for nd in net.hybrid)
    hyb_cys = Dict{Int, Vector{Int}}() ## store the nd numbers of the cycle 
    blobs=biconnectedComponents(net,true)

    desc_mat = descendenceMatrix(net)
    desc_number_dict = Dict(nd.number => findfirst(x-> x==nd.number,desc_mat.nodenumbers_toporder) for nd in net.node) ##map node numbers to their desc ordering
    hyb_number2nd = Dict(nd.number => nd for nd in net.node)

    for (blob,blob_no) in zip(blobs,1:length(blobs))
        nds = unique((x-> getchild(x)).(blob))
        blob_size = length(nds)
        blob_hyb_nds = filter(x->x.hybrid,nds)
        blob_hyb_nds_desc_inds = (x-> findfirst(y-> y==x.number,desc_mat.nodenumbers_toporder)).(blob_hyb_nds)
        blob_level = length(blob_hyb_nds)
        for (hyb,hyb_ind) in zip(blob_hyb_nds,blob_hyb_nds_desc_inds)
            mkpath(output_net_dir*"$(hyb.name)/")
            upper_stack =  filter(x-> ((0 < desc_mat.V[x,hyb_ind]) & (hyb_ind!=x)),blob_hyb_nds_desc_inds) ## all hybs ancestor to hyb_ind
            lower_stack =  filter(x-> ((0 < desc_mat.V[hyb_ind,x]) & (hyb_ind!=x)),blob_hyb_nds_desc_inds) ## all hybs that descend from hyb_ind
            external = isempty(lower_stack)
            stacked = !(isempty(lower_stack)) | !(isempty(upper_stack))

            all_stack_names = (x-> hyb_number2nd[x].name).(desc_mat.nodenumbers_toporder[[upper_stack;lower_stack]]) ## all nodes that share the stack of hyb
            CSV.write(output_net_dir*"$(hyb.name)/stacked.csv",
                DataFrame(stacked_nodes=all_stack_names,
                        islower=[fill(0,length(upper_stack));
                        fill(1,length(lower_stack))])
            )
            

            ##Cycle info for node
            par_nd_numbers = [desc_number_dict[par_nd.number] for par_nd in getparents(hyb)]
            lpar=par_nd_numbers[1]
            rpar=par_nd_numbers[2]
            #Find which nodes have the left and right parents as descendents 
            has_left = [(desc_mat.V[rpar,nd]>0) for nd in 1:length(net.node)]
            has_right= [(desc_mat.V[lpar,nd]>0) for nd in 1:length(net.node)]

            left_cycle =Int[]
            right_cycle=Int[]
            cycle_root=-1
            for (i,(l,r)) in enumerate(zip(reverse(has_left),reverse(has_right))) #go thru has_left and has_right
                i=length(has_right)-i+1
                if (l && r) #IF nodenumbers_toporder is correct, the first time we find a node with both parents as descendents, it should be the 'root' of the cycle
                    cycle_root=i
                    break #we have traversed the cycle if we made it to the 'root'
                elseif l
                    push!(left_cycle,i)
                elseif r
                    push!(right_cycle,i)
                end
            end
            (cycle_root==-1) && error("we didn't find the 'root' of the hybrid cycle starting at node number $(hyb_nd.number)")
            
            cy=[[desc_number_dict[hyb.number]];[cycle_root];left_cycle;right_cycle]
            cy_number = desc_mat.nodenumbers_toporder[cy]
            cy_nds = [hyb_number2nd[k] for k in cy_number]
            cy_hybs = filter(x-> (x.hybrid & (x.number != hyb.number)) ,cy_nds)
            hyb_cys[hyb.number] = cy_number
            for cy_hyb in cy_hybs
                push!(hyb_cycle_map[cy_hyb.number],(hyb.name,true)) #true means that key is a part of value's cycle (upper)
                push!(hyb_cycle_map[hyb.number],(cy_hyb.name,false)) #false means that value is a part of key's  cycle (lower)
            end
            CSV.write(output_net_dir*"$(hyb.name)/cycle_nodes.csv",DataFrame(node_number=cy_number))

            cycle_size = length(cy)
            cycle_diam = sum( (x ->x.length ).(getparentedge.(cy_nds[setdiff(1:end, 2)])))
            cycle_size==2 && (cycle_diam*=2)
            hyb_blob_map[hyb.number]=(blob_size,blob_level,blob_no,external,stacked,cycle_size,cycle_diam)
        end
    end

    for hyb in hyb_nds

        dist50 = abs(0.5-getparentedge(hyb).gamma)
        blob_size,blob_lvl,blob_no,ext,stacked,cycle_size,cycle_diam = hyb_blob_map[hyb.number]

        hyb_cycled = hyb_cycle_map[hyb.number]
        if isempty(hyb_cycled)
            cycled=false
            hyb_cy=DataFrame()
        else
            hyb_cy=DataFrame(hyb_cycle_map[hyb.number],[:hyb,:ispartof])
            cycled=true
        end
        CSV.write(output_net_dir*"/$(hyb.name)/cycle_hybrids.csv",hyb_cy)

        cy_num = hyb_cys[hyb.number]
        cy_adj = false
        for hyb2 in hyb_nds ## check if hybrids share a cycle edge
            (hyb.number == hyb2.number) && continue ## don't check the hyb with itself
            cy_num2 = hyb_cys[hyb2.number]
            if !isempty(intersect(cy_num,cy_num2))
                cy_adj=true
                break #we can stop looking
            end
        end

        push!(hyb_dat,(
            phy_no, 
            hyb.name, #name,
            blob_no,
            dist50, #dist50 - gamma
            blob_size, #blob size
            blob_lvl, #blob level
            ext, #external hyb
            stacked,
            cycle_size,
            cycle_diam,
            cycled,
            cy_adj
        ))
        ##TODO find adjacent cycle hybs - hybrids that aren't a direct part of each other's cycle but share nodes in their cycles.



    end
    




end

bl_dat = DataFrame(bl = branch_lengths)
CSV.write(output_dir*"branch_lengths.csv",bl_dat)
CSV.write(output_dir*"hyb_dat.csv",hyb_dat)
end



art_points = Float64[]
for par_no in 1:18
    output_dir = "../data/pars_$par_no/"
    println("we are at par $par_no" ) 
for phy_no in 1:150
    
    true_net_loc = "../data/pars_$par_no/net_$phy_no/"
    net = readnewick(true_net_loc*"network.extnewick")
    preorder!(net)
    PhyloNetworks.process_biconnectedcomponents!(net)

    for b in net.partition
        if length(b.edges)>=2 #not a trivial bcc
            push!(art_points,PhyloNetworks.number_exitnodes(b)) #number of articulation points
        end        
    end
end
end
art_table = counter(art_points)