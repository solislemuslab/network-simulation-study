using Base.Threads
using PhyloNetworks
using SNaQ
using DataFrames
using CSV

job_id = parseARGS[1]

net = readnewick("network.extnewick")
shrink2cycles!(net)
shrink3cycles!(net)

subnets,subnet_table = PhyloNetworks.getdisplayednetworks(net)
nrows = size(subnet_table)[1]

subnet_levels = repeat([-1],nrows)
for (i,subnet) in enumerate(subnets)
    subnet_levels[i] = getlevel(subnet)
end



liks = repeat([-Inf],nrows)
liks_optimized = repeat([-Inf],nrows)
cfs = readtableCF("CFs.csv")

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
CSV.write("liks.csv",subnet_liks)

