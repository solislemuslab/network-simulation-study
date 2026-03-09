

using CondaPkg
using PythonCall

#we need phylox and physquirrel from Python
const parser = pyimport("phylox.newick_parser")
const psq = pyimport("physquirrel")
const nx = pyimport("networkx")


function newick_to_squirrel_semi_network(net::HybridNetwork)
    
    n = deepcopy(net)
    removedegree2nodes!(n)

    edge_list = Py[]
    ret_edge_list = Py[]
    leaves = String[]
    for e in n.edge
        par = string(getparent(e).number)
        child_nd = getchild(e)
        if child_nd.leaf
            child = child_nd.name
            push!(leaves,child_nd.name)
        else
            child = string(child_nd.number)
        end
        graph_edge = pytuple((par,child))
        push!(edge_list, graph_edge)
        e.hybrid && push!(ret_edge_list,graph_edge)
    end

        network = psq.SemiDirectedNetwork(
        incoming_graph_data = edge_list,  # Pass ALL edges here
        directed_edges = ret_edge_list,     # Pass ONLY directed ones here
        leaves = leaves
    )
    return network
end

function consistency_score(qnets1,qnets2)
    return pyconvert(Float64, qnets1.consistency(qnets2))
end



function hardwiredclusterdistance_firstrooted!(
    root_net::HybridNetwork,
    unroot_net::HybridNetwork,
    taxa::Vector{String};
    M1::Array{Int64, 2})
    net1=deepcopy(root_net)
    net2=deepcopy(unroot_net)

    directedges!(net2) #update containsRoot on edges for possible rootings
    suppressroot!(net2) #remove root unpredictable order
   
    # find all permissible positions for the root
    #=
    net2rootedges = [n for n in net2.edge ]
    for i in length(net2rootedges):-1:1
        println(i)
        try
            rootonedge!(net2, net2rootedges[i])
            suppressroot!(net2)
        catch e
            isa(e, PhyloNetworks.RootMismatch) || rethrow(e)
            deleteat!(net2rootedges, i)
        end
    end
    =#
    net2rootedges = [n for n in net2.edge if n.containroot]
    bestdissimilarity = typemax(Int)

    for n2 in net2rootedges
        #println(n2.number)
        suppressroot!(net2)
        rootonedge!(net2, n2)
        
        ######
        length(setdiff(taxa, String[net2.leaf[i].name for i in 1:net2.numtaxa])) == 0 ||
            error("net1 and net2 do not share the same taxon set. Please prune networks first.")


        ismissing(M1) && (M1 = hardwiredclusters(net1, taxa)) # last row: 10/11 if tree/hybrid edge.
        M2 = hardwiredclusters(net2, taxa)
        dis = 0
        n2ci = collect(1:size(M2, 1)) # cluster indices
        for i1 in 1:size(M1,1)
            found = false
            m1 = 1 .- M1[i1,2:end] # going to the end: i.e. we want to match a tree edge with a tree edge
                                    # and hybrid edge with hybrid edge
            for j in length(n2ci):-1:1 # check only unmatched cluster indices, in reverse
                i2 = n2ci[j]
                if (M1[i1,2:end] == M2[i2,2:end])
                    found = true
                    deleteat!(n2ci, j) # a cluster can be repeated
                    break
                end
            end
            if !found
                dis += 1
            end
        end # (size(M1)[1] - dis) edges have been found in net2, dis edges have not.
        # so size(M2)[1] - (size(M1)[1] - dis) edges in net2 are not in net1.
        diss = dis + dis + size(M2)[1] - size(M1)[1]

        ######



        if diss < bestdissimilarity
            bestdissimilarity = diss
        end
    end
    return bestdissimilarity
end


function only_trees(net::HybridNetwork)
    return net.numhybrids == 0 
end