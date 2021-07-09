#!/usr/bin/julia
#print("hi Julia21rg")
#print("ayuda")
#
using PhyloNetworks, PhyloPlots
using RCall

net0 = readTopology("RNetwork")

R"pdf(file='NetworkjuliaR2.pdf')"
plot(net0, :R)
R"dev.off()"


#savefig("plotJu.png")
print(hybridlambdaformat(net0))
#cd("/home/acosta/Escritorio/Working/Filogenetica/Bash")
#readdir()