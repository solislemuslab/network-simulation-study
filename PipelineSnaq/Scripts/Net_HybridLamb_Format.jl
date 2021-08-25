#print("hi Julia21rg")
#print("ayuda")
#
IniPat = pwd()
#cd(joinpath((IniPat),"RNetworks"))
cd(ARGS[1])
#cd("/home/acosta/GitProjects/network-simulation-study/BashScript_Network_NotLevelOne/RNetworks")
using PhyloNetworks, PhyloPlots
using RCall

#net0 = readTopology("RNetwork")
Nets=split(read(`ls`, String))


for i in 1:length(Nets)
	cd(ARGS[1])
	RNet = Nets[i]
	net0 = readTopology(RNet)
	namefileJu = replace(RNet,"RNetwork"=>"RJuliaNet")*".txt"
	namefileHy = replace(RNet,"RNetwork"=>"JulHybrLamb")
	cd(ARGS[2])
	writeTopology(net0,namefileJu)
	cd(ARGS[3])
	HyLaF = hybridlambdaformat(net0)
	open(namefileHy, "w") do io
		write(io, HyLaF)
	end
	cd(ARGS[1])
end 
