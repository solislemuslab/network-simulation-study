#!/usr/bin/julia
#print("hi Julia21rg")
#print("ayuda")
#
using PhyloNetworks, PhyloPlots
using RCall

#net0 = readTopology("RNetwork")
#cd("/home/acosta/GitProjects/network-simulation-study/BashScript_Network_NotLevelOne")
#writeTopology(net0, "JuliaNetwork1.txt")
#R"pdf(file='NetworkjuliaR2.pdf')"
#plot(net0, :R)
#R"dev.off()"

nets_not_level1 = readMultiTopology("nets_not_level1.extnex")

# bulk format conversion of networks read in nets_not_level1
open("nets_not_level1.hybridlambda", "w") do io    
    for i = 1:length(nets_not_level1)
        hyblambFormatted = hybridlambdaformat(nets_not_level1[i])
        write(io, "$hyblambFormatted\n")
        flush(io)
    end
end
    

#savefig("plotJu.png")
#print(hybridlambdaformat(net0))
#cd("/home/acosta/Escritorio/Working/Filogenetica/Bash")
#readdir()
