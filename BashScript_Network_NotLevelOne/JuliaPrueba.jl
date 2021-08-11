#!/usr/bin/julia
#print("hi Julia21rg")
#print("ayuda")
#
IniPat = pwd()
cd(joinpath((IniPat),"RNetworks"))
#cd("/home/acosta/GitProjects/network-simulation-study/BashScript_Network_NotLevelOne/RNetworks")
using PhyloNetworks, PhyloPlots
using RCall

Nets=split(read(`ls`, String))

#net0 = readTopology("RNetwork")
#cd("/home/acosta/GitProjects/network-simulation-study/BashScript_Network_NotLevelOne")
#writeTopology(net0, "JuliaNetwork1.txt")
#R"pdf(file='NetworkjuliaR2.pdf')"
#plot(net0, :R)
#R"dev.off()"


#savefig("plotJu.png")
#print(hybridlambdaformat(net0))
#cd("/home/acosta/Escritorio/Working/Filogenetica/Bash")
#readdir()

#



for i in 1:length(Nets)
	cd(joinpath((IniPat),"RNetworks"))
	RNet = Nets[i]
	net0 = readTopology(RNet)
	namefileJu = replace(RNet,"RNetwork"=>"RJuliaNet")*".txt"
	namefileHy = replace(RNet,"RNetwork"=>"JulHybrLamb")
	cd(joinpath((IniPat),"RJuliaNetworks"))
	writeTopology(net0,namefileJu)
	cd(joinpath((IniPat),"JuliaForHybLam"))
	HyLaF = hybridlambdaformat(net0)
	open(namefileHy, "w") do io
		write(io, HyLaF)
	end
	cd(joinpath((IniPat),"RNetworks"))
end